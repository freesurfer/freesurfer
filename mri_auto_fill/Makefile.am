## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_auto_fill
mri_auto_fill_SOURCES=mri_auto_fill.c
mri_auto_fill_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_auto_fill_LDFLAGS=$(OS_LDFLAGS)

