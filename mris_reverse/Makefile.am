## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mris_reverse
mris_reverse_SOURCES=mris_reverse.c
mris_reverse_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_reverse_LDFLAGS=$(OS_LDFLAGS)

