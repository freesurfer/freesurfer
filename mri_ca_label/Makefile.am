## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_ca_label
mri_ca_label_SOURCES=mri_ca_label.c
mri_ca_label_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_ca_label_LDFLAGS=$(OS_LDFLAGS)

