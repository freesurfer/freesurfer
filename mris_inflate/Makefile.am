## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mris_inflate
mris_inflate_SOURCES=mris_inflate.c
mris_inflate_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_inflate_LDFLAGS=$(OS_LDFLAGS)

