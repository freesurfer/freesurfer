## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_surfacemask
mri_surfacemask_SOURCES=mri_surfacemask.cpp
mri_surfacemask_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_surfacemask_LDFLAGS=$(OS_LDFLAGS)

