## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_watershed
mri_watershed_SOURCES=mri_watershed.cpp
mri_watershed_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_watershed_LDFLAGS=$(OS_LDFLAGS)

SUBDIRS=brain_volume
