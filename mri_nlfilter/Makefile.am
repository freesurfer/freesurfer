## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_nlfilter
mri_nlfilter_SOURCES=mri_nlfilter.c
mri_nlfilter_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_nlfilter_LDFLAGS=$(OS_LDFLAGS)

