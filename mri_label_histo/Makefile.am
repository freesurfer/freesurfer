## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_label_histo
mri_label_histo_SOURCES=mri_label_histo.c
mri_label_histo_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_label_histo_LDFLAGS=$(OS_LDFLAGS)

