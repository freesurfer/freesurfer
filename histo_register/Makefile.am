# 
# Makefile.am
#

if BUILDOPENCV
SBL=$(top_srcdir)/histo_register/SimpleBaseLib
PVL=$(top_srcdir)/histo_register/ParticleVideoLib
SUBDIRS=./SimpleBaseLib ./ParticleVideoLib

bin_PROGRAMS = histo_register
histo_register_SOURCES=\
	prep/BlockPrep.cc \
	prep/HistoPrep.cc \
	prep/HistoStats.cc \
	prep/HistoStitch.cc \
	prep/MPrep.cc \
	prep/Polarization.cc \
	prep/Prediction.cc \
	prep/StitchNode.cc \
	prep/VolumeFile.cc \
	prep/VolumeUtil.cc \
	registration/BlockRegister.cc \
	registration/CorresField3D.cc \
	registration/HistoRegister.cc \
	registration/HistoTransform.cc \
	registration/ImageSetSeq.cc \
	registration/Normalization.cc \
	registration/TestCorres3D.cc \
	registration/VarCorres3D.cc \
	registration/VarCorres3DUtil.cc \
	Main.cc

histo_register_CXXFLAGS=-I. -I$(SBL)/include -I$(SBL)/external -I$(PVL)/include $(OPENCV_CXXFLAGS)
histo_register_CPPFLAGS=-DUSE_OPENCV -DUSE_CDT
LOCALLIBS=-L$(SBL) -L$(PVL) -lpvl -lsbl
histo_register_LDADD=$(OPENCV_LIBS) $(LOCALLIBS)


install-exec-hook:
	cp $(top_builddir)/histo_register/histo_register \
		$(DESTDIR)$(bindir)/histo_register.bin
	echo "#!/bin/tcsh -ef" > $(DESTDIR)$(bindir)/histo_register
	echo "setenv LD_LIBRARY_PATH /usr/pubsw/packages/opencv/current/lib" \
		>> $(DESTDIR)$(bindir)/histo_register
	echo "histo_register.bin \$$argv" >> $(DESTDIR)$(bindir)/histo_register
	chmod a+x $(DESTDIR)$(bindir)/histo_register

uninstall-hook:
	rm -f $(DESTDIR)$(bindir)/histo_register.bin

endif
