## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mris_expand
mris_expand_SOURCES=mris_expand.c
mris_expand_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_expand_LDFLAGS=$(OS_LDFLAGS)

