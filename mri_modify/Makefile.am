## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_modify
mri_modify_SOURCES=mri_modify.cpp
mri_modify_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_modify_LDFLAGS=$(OS_LDFLAGS)

