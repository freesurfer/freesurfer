## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mris_rescale
mris_rescale_SOURCES=mris_rescale.c
mris_rescale_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_rescale_LDFLAGS=$(OS_LDFLAGS)

