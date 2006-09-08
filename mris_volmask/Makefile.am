##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom \
$(VTK_CXXFLAGS)/vtk-5.0 -Wno-deprecated

bin_PROGRAMS = mris_volmask
mris_volmask_SOURCES=mris_volmask.cpp cmd_line_interface.cpp
mris_volmask_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH)) \
$(X_LIBS) $(GL_LIBS) \
$(VTK_LIBS) \
-lvtkFiltering -lvtkHybrid -lvtkRendering -lvtkGraphics \
-lvtkImaging -lvtkftgl -lvtkfreetype \
-lvtkIO -lvtkFiltering -lvtkCommon -lvtkDICOMParser \
-lvtkpng -lvtktiff -lvtkzlib -lvtkjpeg -lvtkexpat -lvtksys \
-ldl -lvtkMPEG2Encode -lvtkexoIIc -lvtkNetCDF

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
