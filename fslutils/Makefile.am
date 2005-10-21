## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=

bin_PROGRAMS = fsl_label2voxel
fsl_label2voxel_SOURCES=fsl_label2voxel.cpp
fsl_label2voxel_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
fsl_label2voxel_LDFLAGS=$(OS_LDFLAGS)



# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
