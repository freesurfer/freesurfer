## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_cor2label
mri_cor2label_SOURCES=mri_cor2label.c
mri_cor2label_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_cor2label_LDFLAGS=$(OS_LDFLAGS)

