## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mris_segment_vals
mris_segment_vals_SOURCES=mris_segment_vals.c
mris_segment_vals_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_segment_vals_LDFLAGS=$(OS_LDFLAGS)

