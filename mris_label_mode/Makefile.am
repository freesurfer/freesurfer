## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mris_label_mode
mris_label_mode_SOURCES=mris_label_mode.c
mris_label_mode_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_label_mode_LDFLAGS=$(OS_LDFLAGS)

