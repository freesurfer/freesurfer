## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mris_thickness_comparison
mris_thickness_comparison_SOURCES=mris_thickness_comparison.c
mris_thickness_comparison_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_thickness_comparison_LDFLAGS=$(OS_LDFLAGS)

