## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_synthesize
mri_synthesize_SOURCES=mri_synthesize.c
mri_synthesize_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_synthesize_LDFLAGS=$(OS_LDFLAGS)

