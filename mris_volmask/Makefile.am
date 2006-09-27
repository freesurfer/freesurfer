##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom \
$(VTK_CXXFLAGS)/vtk-5.0 -Wno-deprecated

if ENABLE_VTK_APPS
bin_PROGRAMS = mris_volmask
EXTRA_DIST = cmd_line_interface.h
else
bin_PROGRAMS =
EXTRA_DIST = 
endif
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

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
