## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_parselabel
mri_parselabel_SOURCES=mri_parselabel.cpp
mri_parselabel_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_parselabel_LDFLAGS=$(OS_LDFLAGS)

