## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_normalize
mri_normalize_SOURCES=mri_normalize.c
mri_normalize_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_normalize_LDFLAGS=$(OS_LDFLAGS)

