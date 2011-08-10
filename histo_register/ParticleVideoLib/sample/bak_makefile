CPP=g++
OPTLEVEL=2
WARNFLAGS=-Wall -Wstrict-aliasing -Werror
INCLUDE=-I. -I$(SBL)/include -I$(PVL)/include -I$(OPENCV)/include
OPTFLAGS=-msse2 -mfpmath=sse -march=i686
CONFIGFLAGS=-DUSE_OPENCV -DUSE_CDT #-D_DEBUG 
CPPFLAGS=-g -fpack-struct=8 -mpc64 -fno-strict-aliasing -ffloat-store 
CPPFLAGS+=$(WARNFLAGS) $(INCLUDE) -O$(OPTLEVEL) $(OPTFLAGS) $(CONFIGFLAGS) # -pg 
SOURCES=main.cc
OBJECTS=$(SOURCES:.cc=.o)
EXECUTABLE=pvsample
LIBPATH=-L$(SBL)/lib -L$(PVL)/lib
LIBRARIES=-lm -lopencv_core -lopencv_imgproc -lopencv_highgui -lopencv_ml -lpvl -lsbl

.PHONY: depend clean

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CPP) $(CPPFLAGS) $(OBJECTS) -o $(EXECUTABLE) $(LIBPATH) $(LIBRARIES)

.cc:
	$(CPP) $(CPPFLAGS) -c $@.cc -o $@.o

clean:
	rm -f *.o $(EXECUTABLE)

depend: $(SOURCES)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it
