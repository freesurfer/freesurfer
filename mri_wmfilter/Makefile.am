## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_wmfilter
mri_wmfilter_SOURCES=mri_wmfilter.c
mri_wmfilter_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_wmfilter_LDFLAGS=$(OS_LDFLAGS)

