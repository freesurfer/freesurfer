## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_extract
mri_extract_SOURCES=mri_extract.c
mri_extract_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_extract_LDFLAGS=$(OS_LDFLAGS)

