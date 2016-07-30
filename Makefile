CFITSIO = $(FITSIOROOT)
CPP = g++
CC = gcc
CFLAGS = -Wall -I$(CFITSIO) $(shell root-config --cflags)
LIBS = -L$(CFITSIO) -lcfitsio -lm $(shell root-config --libs)
GLIBS = 
GLIBS += 
OBJECTS = skipper2root.o 
HEADERS = globalConstants.h

ALL : skipper2root.exe
	@echo "Listo!"

skipper2root.exe : $(OBJECTS)
	$(CPP) $(OBJECTS) -o skipper2root.exe $(LIBS) $(GLIBS) $(CFLAGS)

skipper2root.o : skipper2root.cc $(HEADERS)
	$(CPP) -c skipper2root.cc -o skipper2root.o $(CFLAGS)

clean:
	rm -f *~ *.o *.exe