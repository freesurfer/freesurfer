## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mris_classify_thickness
mris_classify_thickness_SOURCES=mris_classify_thickness.c
mris_classify_thickness_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_classify_thickness_LDFLAGS=$(OS_LDFLAGS)

