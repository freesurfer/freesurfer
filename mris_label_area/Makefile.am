## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mris_label_area
mris_label_area_SOURCES=mris_label_area.c
mris_label_area_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_label_area_LDFLAGS=$(OS_LDFLAGS)

