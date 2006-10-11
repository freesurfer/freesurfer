##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include $(VTK_CXXFLAGS) -Wno-deprecated

if ENABLE_VTK_APPS
bin_PROGRAMS = mris_volmask
EXTRA_DIST = cmd_line_interface.h

mris_volmask_SOURCES=mris_volmask.cpp cmd_line_interface.cpp
mris_volmask_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH)) \
$(X_LIBS) $(GL_LIBS) \
$(VTK_LIBS) \
-ldl -lpthread \
-lvtkFiltering -lvtkHybrid -lvtkRendering -lvtkGraphics \
-lvtkImaging -lvtkftgl -lvtkfreetype \
-lvtkIO -lvtkFiltering -lvtkCommon -lvtkDICOMParser \
-lvtkpng -lvtktiff -lvtkzlib -lvtkjpeg -lvtkexpat -lvtksys \
-lvtkMPEG2Encode -lvtkexoIIc -lvtkNetCDF

# put a wrapper around the bin, used to setup VTK enviro vars
install-exec-hook:
	mv $(bindir)/mris_volmask $(bindir)/mris_volmask.bin
	echo "#!/bin/tcsh -ef" > $(bindir)/mris_volmask
	echo "source \$$FREESURFER_HOME/bin/vtk_setup" >> $(bindir)/mris_volmask
	echo "mris_volmask.bin \$$argv" >> $(bindir)/mris_volmask
	chmod a+x $(bindir)/mris_volmask

else
bin_PROGRAMS =
EXTRA_DIST = 
endif

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
