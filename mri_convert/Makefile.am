## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard
AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_convert
mri_convert_SOURCES=mri_convert.c
mri_convert_LDADD= $(addprefix $(top_srcdir)/,$(LIBS_MGH))
mri_convert_DEPENDENCIES= $(addprefix $(top_srcdir)/,$(LIBS_MGH))

