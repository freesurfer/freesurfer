## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mris_compute_overlap
mris_compute_overlap_SOURCES=mris_compute_overlap.c
mris_compute_overlap_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_compute_overlap_LDFLAGS=$(OS_LDFLAGS)

