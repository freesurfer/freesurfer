## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mris_convert
mris_convert_SOURCES=mris_convert.c
mris_convert_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_convert_LDFLAGS=$(OS_LDFLAGS)

