## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=

bin_PROGRAMS = mri_annotation2label
mri_annotation2label_SOURCES=mri_annotation2label.c
mri_annotation2label_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_annotation2label_LDFLAGS=$(OS_LDFLAGS)

