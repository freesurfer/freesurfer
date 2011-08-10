SBL=../../SimpleBaseLib-0.0.4
PVL=../../ParticleVideoLib-0.0.2
OPENCV=/usr/pubsw/packages/opencv/current
CPP=g++
OPTLEVEL=2
WARNFLAGS=-Wall -Wstrict-aliasing -Werror
INCLUDE=-I../include -I$(OPENCV)/include -I$(SBL)/include
OPTFLAGS=-msse2 -mfpmath=sse
CONFIGFLAGS=-DUSE_OPENCV -DUSE_CDT #-D_DEBUG 
CPPFLAGS=-fpack-struct=8
CPPFLAGS+=$(WARNFLAGS) $(INCLUDE) -O$(OPTLEVEL) $(OPTFLAGS) $(CONFIGFLAGS) # -g -pg 
SOURCES=$(wildcard ../src/*.cc) 
INCLUDES=$(wildcard ../include/pvl/*.h)
OBJECTS=$(SOURCES:.cc=.o)
LIBRARY=../lib/libpvl.a

.PHONY: depend clean valgrind

all: $(SOURCES) $(LIBRARY)

$(LIBRARY): $(OBJECTS)
	ar rcs $(LIBRARY) $(OBJECTS)

.cc:
	$(CPP) $(CPPFLAGS) -c $@.cc -o $@.o

clean:
	rm -f ../src/*.o $(LIBRARY)

depend: $(SOURCES)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it
