##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include $(VTK_CXXFLAGS) -Wno-deprecated

bin_PROGRAMS = mris_volmask

mris_volmask_SOURCES=mris_volmask.cpp cmd_line_interface.cpp
mris_volmask_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH)) \
$(X_LIBS) $(GL_LIBS) \
-ldl -lpthread \
$(VTK_LIBS) \
$(VTK_FREETYPE_LIB) $(VTK_PNG_LIB) \
$(VTK_TIFF_LIB) $(VTK_JPEG_LIB) \
$(VTK_Z_LIB) $(VTK_EXPAT_LIB) \
-lvtkFiltering $(VTK_VERDICT_LIB) -lvtkHybrid \
-lvtkRendering -lvtkGraphics \
-lvtkImaging -lvtkftgl \
-lvtkIO -lvtkFiltering -lvtkCommon -lvtkDICOMParser \
-lvtksys \
$(VTK_MPEG2ENCODE_LIB) -lvtkexoIIc -lvtkNetCDF

# put a wrapper around the bin, used to setup VTK enviro vars
install-exec-hook:
	cp $(srcdir)/mris_volmask $(DESTDIR)$(bindir)/mris_volmask.bin
	echo "#!/bin/tcsh -ef" > $(DESTDIR)$(bindir)/mris_volmask
	echo "source \$$FREESURFER_HOME/bin/vtk_setup" \
		>> $(DESTDIR)$(bindir)/mris_volmask
	echo "mris_volmask.bin \$$argv" >> $(DESTDIR)$(bindir)/mris_volmask
	chmod a+x $(DESTDIR)$(bindir)/mris_volmask

EXTRA_DIST=cmd_line_interface.h

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
