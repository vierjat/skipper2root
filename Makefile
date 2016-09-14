CFITSIO = $(FITSIOROOT)
CPP = g++
CC = gcc
GCCNEWERTHAN47 := $(shell expr `gcc -dumpversion` \>= 4.7)
CFLAGS = -Wall -I$(CFITSIO) $(shell root-config --cflags)
LIBS = -L$(CFITSIO) -lcfitsio -lm $(shell root-config --libs)
GLIBS = 
GLIBS +=
ifeq "$(GCCNEWERTHAN47)" "1"
  CFLAGS += -std=c++11
else
  CFLAGS += -std=c++0x
endif
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
