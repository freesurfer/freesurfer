## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mris_flatten
mris_flatten_SOURCES=mris_flatten.c
mris_flatten_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_flatten_LDFLAGS=$(OS_LDFLAGS)

