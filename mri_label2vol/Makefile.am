## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_label2vol
mri_label2vol_SOURCES=mri_label2vol.c
mri_label2vol_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_label2vol_LDFLAGS=$(OS_LDFLAGS)

