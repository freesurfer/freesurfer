## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_fill
mri_fill_SOURCES=mri_fill.c
mri_fill_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_fill_LDFLAGS=$(OS_LDFLAGS)

