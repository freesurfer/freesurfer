## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_brain_volume
mri_brain_volume_SOURCES=mri_brain_volume.cpp
mri_brain_volume_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_brain_volume_LDFLAGS=$(OS_LDFLAGS)

