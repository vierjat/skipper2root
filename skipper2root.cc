#include <string.h>
#include <stdio.h>
#include "fitsio.h"

#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cerrno>
#include <memory>

#include <sys/time.h>
#include <time.h>
#include <inttypes.h>
#include <fstream>
#include <unistd.h>
#include <getopt.h>    /* for getopt_long; standard getopt is in unistd.h */
#include <vector>
#include <algorithm>
#include <ctime>
#include <climits>
#include <cmath>
#include <iomanip>
#include <numeric>

#include "TFile.h"
#include "TTree.h"

#include "globalConstants.h"

using namespace std;

const long imCols    = 1000;
const long ccdCols   = 452*2;
const long nIgnore   = 3;
const long nOS       = (imCols/2 - ccdCols/2)-nIgnore-1;
const long lOsStart  = ccdCols/2+nIgnore; // Left OS start
const long lOsEnd    = lOsStart+nOS;        // Left OS end
const long rOsStart  = imCols/2;          // Right OS start
const long rOsEnd    = rOsStart+nOS;        // Right OS end

const int  nMeanTrim = 4;

const int kSaveSamplesFlag = 1;
const int kZeroThrFlag     = 2;

const int kMaxNSpl = 10000;

int deleteFile(const char *fileName){
  cout << yellow;
  cout << "Will overwrite: " << fileName << endl << endl;
  cout << normal;
  return unlink(fileName);
}

bool fileExist(const char *fileName){
  ifstream in(fileName,ios::in);
  
  if(in.fail()){
    //cout <<"\nError reading file: " << fileName <<"\nThe file doesn't exist!\n\n";
    in.close();
    return false;
  }
  
  in.close();
  return true;
}

/*========================================================
  ASCII progress bar
==========================================================*/
void showProgress(unsigned int currEvent, unsigned int nEvent) {

  const int nProgWidth=50;

  if ( currEvent != 0 ) {
    for ( int i=0;i<nProgWidth+8;i++)
      cout << "\b";
  }

  double percent = (double) currEvent/ (double) nEvent;
  int nBars = (int) ( percent*nProgWidth );

  cout << " |";
  for ( int i=0;i<nBars-1;i++)
    cout << "=";
  if ( nBars>0 )
    cout << ">";
  for ( int i=nBars;i<nProgWidth;i++)
    cout << " ";
  cout << "| " << setw(3) << (int) (percent*100.) << "%";
  cout << flush;

}

void printCopyHelp(const char *exeName, bool printFullHelp=false){
  
  if(printFullHelp){
    cout << bold;
    cout << endl;
    cout << "This program process the raw Skipper CCD data. It computes overscan\n"
         << "mean for each sample and subtracts it line by line.\n"
         << "The output file will be a ROOT image containing a TTree with the pixels\n"
         << "values averaged over all the samples (after subtraction of the corresponding\n"
         << "overscan value). It's also possible to save an additional TTree that will\n"
         << "contain the individual values of all the samples.\n";
    cout << normal;
  }
  cout << "==========================================================================\n";
  cout << red;
  cout << "\nUsage:\n";
  cout << "  "   << exeName << " <input file> -o <output filename> \n\n";
  cout << "\nOptions:\n";
  cout << "  -s for saving the individual values of all the samples.\n";
  cout << "  -d for overwriting the output file if it exist.\n\n";
  cout << "  -z <zero threshold in ADC> for using pixels with skPix vale smaller than zeroThr in the OS mean.\n\n";
  cout << normal;
  cout << blue;
  cout << "For any problems or bugs contact Javier Tiffenberg <javiert@fnal.gov>\n\n";
  cout << normal;
  cout << "==========================================================================\n\n";
}

int procSkipperImage(const char *inFile, const char *outF, const int opt = 0, const double zeroThr = 100000){

  const bool saveSamples = (opt & kSaveSamplesFlag);
  const bool useZeroThr  = (opt & kZeroThrFlag);

  /* Do not overwrite the output file if it already exist */
  if(fileExist(outF)){
    cout << yellow << "\nThe output file exist.\n" << normal;
    cout << red    << "Will NOT continue.\n" << normal;
    return -10;
  }

  fitsfile *infptr;
  
  int status = 0;  
  int nkeys;
  int nhdu = 0;
  long totpix = 0;
  //char card[81];
  
  fits_open_file(&infptr, inFile, READONLY, &status); /* Open the input file */
  if (status != 0) return(status);
  
  fits_get_num_hdus(infptr, &nhdu, &status); // get the number of HDUs
  
  TFile outRootFile(outF,"RECREATE");


  Int_t x;
  Int_t y;
  Int_t ohdu;
  Double_t pix;
  Double_t osMean;

  TTree skPixTree("skPixTree", "skPixTree");
  skPixTree.Branch("x",    &x, "x/I");
  skPixTree.Branch("y",    &y, "y/I");
  skPixTree.Branch("ohdu", &ohdu, "ohdu/I");
  skPixTree.Branch("pix",  &pix, "pix/D");
  skPixTree.Branch("osMean",  &osMean, "osMean/D");

  Int_t nSpl;
  Double_t splPix[kMaxNSpl];
  TTree splPixTree("splPixTree", "splPixTree");
  splPixTree.Branch("x",    &x, "x/I");
  splPixTree.Branch("y",    &y, "y/I");
  splPixTree.Branch("ohdu", &ohdu, "ohdu/I");
  splPixTree.Branch("nSpl",  &nSpl, "nSpl/I");
  splPixTree.Branch("pix",   &(splPix[0]), "pix[nSpl]/D");
  splPixTree.Branch("skPix", &pix, "skPix/D");

  for(int eI=1; eI<=nhdu; ++eI){  /* Main loop through each extension */
    
    int hdutype, bitpix, naxis = 0;
    long naxes[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
    fits_movabs_hdu(infptr, eI, &hdutype, &status);
    
    fits_get_img_param(infptr, 9, &bitpix, &naxis, naxes, &status); /* get image dimensions and total number of pixels in image */
    totpix = naxes[0] * naxes[1] * naxes[2] * naxes[3] * naxes[4] * naxes[5] * naxes[6] * naxes[7] * naxes[8];
    
    if (hdutype != IMAGE_HDU || naxis == 0 || totpix == 0){
      /* ignore */
    }
    else{
      fits_get_hdrspace(infptr, &nkeys, NULL, &status);
      long fpixel[2] = {1,1};
      double nulval = 0.;
      int anynul = 0;
      double* inArray = new double[totpix];
      fits_read_pix(infptr, TDOUBLE, fpixel, totpix, &nulval, inArray, &anynul, &status);


      /* create output image */
      const int nSamples = naxes[0]/imCols;
      if(nSamples > kMaxNSpl){
        cerr << red << "\nERROR: nSamples > kMaxNSpl !!!\n\nWill not continue.\n\n;" <<normal << endl;
        break;
      }
      long totpixSkp = totpix/nSamples;
      double* outArray = new double[totpixSkp];
      double* outMeanArray = new double[totpixSkp];
      Double_t* fullOutArray = new Double_t[totpix];
      for (int j = 0; j < totpixSkp; ++j) outMeanArray[j] = 0;
      naxes[0] /= nSamples;
      const int nRows = naxes[1];
      ohdu = eI;
      

      /* These are only used if kZeroThrFlag is set */ 
      vector< vector<double> > lSkOsFirstPass(nRows, std::vector<double>(nOS,0) );
      vector< vector<double> > rSkOsFirstPass(nRows, std::vector<double>(nOS,0) );

      if(useZeroThr){
        for (int s = 0; s < nSamples; ++s){
          double lOsSplV[nOS];
          double rOsSplV[nOS];
          for (int r = 0; r < nRows; ++r){
            for (long o = 0; o < nOS; ++o){
              lOsSplV[o]=inArray[nSamples*(o+r*imCols+lOsStart)+s];
              rOsSplV[o]=inArray[nSamples*(o+r*imCols+rOsStart)+s];
            }

            /* compute stable mean for the OS pixels */
            double osV[nOS];
            // left side
            partial_sort_copy(lOsSplV, lOsSplV+nOS, osV, osV+nOS);
            double lMean = (nOS>2*nMeanTrim)? accumulate(osV+nMeanTrim, osV+(nOS-nMeanTrim), 0.0) : accumulate(osV, osV, 0.0);
            lMean = (nOS>4)? lMean/(nOS-nMeanTrim*2) : lMean/nOS;
            // right side
            partial_sort_copy(rOsSplV, rOsSplV+nOS, osV, osV+nOS);
            double rMean = (nOS>2*nMeanTrim)? accumulate(osV+nMeanTrim, osV+(nOS-nMeanTrim), 0.0) : accumulate(osV, osV, 0.0);
            rMean = (nOS>4)? rMean/(nOS-nMeanTrim*2) : rMean/nOS;
            
            // subtract OS mean for each sample in the OS
            for (int o = 0; o < nOS; ++o){
              lSkOsFirstPass[r][o] += lOsSplV[o] - lMean;
              rSkOsFirstPass[r][o] += rOsSplV[o] - rMean;
            }
          }
        }
     
        for (int r = 0; r < nRows; ++r){
          for (long o = 0; o < nOS; ++o){
            lSkOsFirstPass[r][o] /= nSamples;
            rSkOsFirstPass[r][o] /= nSamples;
          }
        }
      }
      /* End of first look in zeroThr processing */


      /* Subtract OS mean for each sample */
      double osAuxV[nOS];
      for (int s = 0; s < nSamples; ++s){ // loop on samples

        for (long j = 0; j < totpixSkp; ++j) outArray[j]=inArray[nSamples*j+s];

        for (int r = 0; r < nRows; ++r){ // loop on rows
          double* rowPtr = outArray+imCols*r;
          double lMean = 0;
          double rMean = 0;

          if(useZeroThr){
            int lNZero = 0;
            int rNZero = 0;
            for (long o = 0; o < nOS; ++o){
              if(lSkOsFirstPass[r][o]<zeroThr){
                lMean += rowPtr[lOsStart+o];
                ++lNZero;
              }
              if(rSkOsFirstPass[r][o]<zeroThr){
                rMean += rowPtr[rOsStart+o];
                ++rNZero;
              }
            }

            lMean /= lNZero;
            rMean /= rNZero;
          }
          else{
            /* compute stable mean for the OS pixels */
            // left side
            partial_sort_copy(rowPtr+lOsStart, rowPtr+lOsEnd, osAuxV, osAuxV+nOS);
            lMean = (nOS>2*nMeanTrim)? accumulate(osAuxV+nMeanTrim, osAuxV+(nOS-nMeanTrim), 0.0) : accumulate(osAuxV, osAuxV, 0.0);
            lMean = (nOS>4)? lMean/(nOS-nMeanTrim*2) : lMean/nOS;
            // right side
            partial_sort_copy(rowPtr+rOsStart, rowPtr+rOsEnd, osAuxV, osAuxV+nOS);
            rMean = (nOS>2*nMeanTrim)? accumulate(osAuxV+nMeanTrim, osAuxV+(nOS-nMeanTrim), 0.0) : accumulate(osAuxV, osAuxV, 0.0);
            rMean = (nOS>4)? rMean/(nOS-nMeanTrim*2) : rMean/nOS;
          }

          // subtract OS mean for each pixel on the row for the current sample
          for (int c = 0; c < imCols/2; ++c){ // loop on the pixels on the row
            *(rowPtr+c) -= lMean;
            *(rowPtr+imCols/2+c) -= rMean;
          }

        }
        for (long j = 0; j < totpixSkp; ++j) outMeanArray[j] += outArray[j];
        for (long j = 0; j < totpixSkp; ++j) fullOutArray[j*nSamples+s] = outArray[j];

        if(gVerbosity) showProgress(s + (eI-1)*nSamples, nhdu*nSamples);
      }

      if(saveSamples){
        nSpl = nSamples;
        for (long j = 0; j < totpixSkp; ++j){
          x = j%imCols;
          y = j/imCols;
          pix = 0;
          for (int s = 0; s < nSamples; ++s){
            splPix[s] = fullOutArray[j*nSamples + s];
            pix += splPix[s];
          }
          pix /= nSamples;
          splPixTree.Fill();
        }
        splPixTree.Write();
      }

      for (int j = 0; j < totpixSkp; ++j) outMeanArray[j] /= nSamples;

      for (long j = 0; j < totpixSkp; ++j){
        x = j%imCols;
        y = j/imCols;
        pix = outMeanArray[j];
        skPixTree.Fill();
      }

      delete[] inArray;
      delete[] outArray;
      delete[] outMeanArray;
      delete[] fullOutArray;
    }
  }

  skPixTree.Write();
  outRootFile.Close();

  fits_close_file(infptr,   &status);
  if(gVerbosity){
    showProgress(1, 1);
  }
  return status;
}


double mean(const double *v, const int &N){
  
  std::vector<double> temparray(v, v+N);
  std::sort(temparray.begin(), temparray.end());
  
  const int nMin = N/3;
  const int nMax = 2*N/3;
  
  double sum=0;
  for(int i=nMin;i<nMax;++i){
    sum+=temparray[i];
  }
  return sum/(nMax-nMin);
}


bool isSaturated(const double &pixVal, const int &bitpix, const double &bzero){
  
  float saturationVal = 0;
  
  switch(bitpix) {
      case BYTE_IMG:
          saturationVal = 128+bzero;
          break;
      case SHORT_IMG:
          saturationVal = 32768+bzero;
          break;
      default:
          saturationVal = kSatValue;
  }
  
  if(pixVal>=saturationVal*kSatMargin)
    return true;
  
  return false;
}

void checkArch(){
  if(sizeof(float)*CHAR_BIT!=32 || sizeof(double)*CHAR_BIT!=64){
    cout << red;
    cout << "\n ========================================================================================\n";
    cout << "   WARNING: the size of the float and double variables is non-standard in this computer.\n";
    cout << "   The program may malfunction or produce incorrect results\n";
    cout << " ========================================================================================\n";
    cout << normal;
  }
}

int processCommandLineArgs(const int argc, char *argv[], string &inFile, string &outFile, int &flags, double &zeroThr){
  
  if(argc == 1) return 1;
  
  zeroThr            = 100000;
  bool outFileFlag   = false;
  bool owOutFileFlag = false;
  int opt=0;
  while ( (opt = getopt(argc, argv, "i:o:z:dsqQhH?")) != -1) {
    switch (opt) {
    case 'o':
      if(!outFileFlag){
        outFile = optarg;
        outFileFlag = true;
      }
      else{
        cerr << red << "\nError, can not set more than one output file!\n\n" << normal;
        return 2;
      }
      break;
    case 'd': /* overwrite output fie is exist */
      owOutFileFlag = true;
      break;
    case 's':
      flags |= kSaveSamplesFlag;
      break;
    case 'z':
      if(zeroThr==100000){
        flags |= kZeroThrFlag;
        char * e;
        errno = 0;
        zeroThr = std::strtod(optarg, &e);
        if (*e != '\0' /* didn't consume the entire string */ || errno != 0 /* error, overflow or underflow */ ){
          cerr << red << "\nError reading zeroThr value.\nWill NOT continue.\n\n"  << normal;
          return 3;
        }
      }
      else{
        cerr << red << "\nError, can not set more than one zeroThr!\nWill NOT continue.\n\n"  << normal;
        return 4;
      }
      break;
    case 'Q':
    case 'q':
      gVerbosity = 0;
      break;
    case 'h':
    case 'H':
    default: /* '?' */
      return 1;
    }
  }
  
  if(!outFileFlag){
    cerr << red << "\nError: output filename missing.\n" << normal;
    return 2;
  }

  inFile="";
  
  if(argc-optind==0){
    cerr << red << "Error: no input file provided!\n\n" << normal;
    return 1;
  }
  else if(argc-optind>1){
    cerr << red << "Error: more than one input file provided!\n\n" << normal;
    return 1;
  }
  
  inFile=argv[optind];
  if(!fileExist(inFile.c_str())){
    cout << red << "\nError reading input file: " << inFile <<"\nThe file doesn't exist!\n\n" << normal;
    return 1;
  }

  if(outFile == inFile){
    cerr << red << "\nError: The output file can't be the input file!\n" << normal;
    cerr << red << "Will NOT continue\n\n" << normal;
    return -102;
  }

  /* Overwrite the output file if it already exist */
  if(fileExist(outFile.c_str())){
    cout << yellow << "\nThe output file exist.\n" << normal;
    if(owOutFileFlag){
      cout << yellow << "Will overwrite the output file.\n\n" << normal;
      deleteFile(outFile.c_str());
    }
    else{
      cout << yellow << "Please provide a different name or use the \"-d\" option.\n\n" << normal;
      return 5;
    }
  }
  
  return 0;
}

int main(int argc, char *argv[])
{
  time_t start,end;
  double dif;
  time (&start);


  string outFile;
  string inFile;
  int opt = 0;
  double zeroThr = -1;
  
  int returnCode = processCommandLineArgs( argc, argv, inFile, outFile, opt, zeroThr);
  if(returnCode!=0){
    if(returnCode == 1) printCopyHelp(argv[0],true);
    if(returnCode == 2) printCopyHelp(argv[0]);
    return returnCode;
  }
  
  if(gVerbosity){
    cout << bold << "\nWill read the following file:\n" << normal;
    cout << "\t" << inFile << endl;
    cout << bold << "\nThe output will be saved in the file:\n\t" << normal << outFile << endl;
  }

  int status  = procSkipperImage(inFile.c_str(), outFile.c_str(), opt, zeroThr);

  if (status != 0){ 
    fits_report_error(stderr, status);
    return status;
  }
  
  /* Report */
  time (&end);
  dif = difftime (end,start);
  if(gVerbosity) cout << green << "\nAll done!\n" << bold << "-> It took me " << dif << " seconds to do it!\n\n" << normal;

  // return status;
  return 0;
}


