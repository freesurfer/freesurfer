## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mris_ca_label
mris_ca_label_SOURCES=mris_ca_label.c
mris_ca_label_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_ca_label_LDFLAGS=$(OS_LDFLAGS)

