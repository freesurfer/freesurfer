## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mris_ms_refine
mris_ms_refine_SOURCES=mris_ms_refine.c
mris_ms_refine_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_ms_refine_LDFLAGS=$(OS_LDFLAGS)

