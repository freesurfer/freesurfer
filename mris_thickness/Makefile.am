## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mris_thickness
mris_thickness_SOURCES=mris_thickness.c
mris_thickness_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_thickness_LDFLAGS=$(OS_LDFLAGS)

