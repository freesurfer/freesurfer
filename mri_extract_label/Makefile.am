## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_extract_label
mri_extract_label_SOURCES=mri_extract_label.c
mri_extract_label_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_extract_label_LDFLAGS=$(OS_LDFLAGS)

