## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mris_sample_label
mris_sample_label_SOURCES=mris_sample_label.c
mris_sample_label_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_sample_label_LDFLAGS=$(OS_LDFLAGS)

