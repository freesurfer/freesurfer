## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = fsl_label2voxel
fsl_label2voxel_SOURCES=fsl_label2voxel.cpp
fsl_label2voxel_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
fsl_label2voxel_LDFLAGS=$(OS_LDFLAGS)

