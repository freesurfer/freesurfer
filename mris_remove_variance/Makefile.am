## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mris_remove_variance
mris_remove_variance_SOURCES=mris_remove_variance.c
mris_remove_variance_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_remove_variance_LDFLAGS=$(OS_LDFLAGS)

