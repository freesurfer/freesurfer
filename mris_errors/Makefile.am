## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mris_errors
mris_errors_SOURCES=mris_errors.c
mris_errors_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_errors_LDFLAGS=$(OS_LDFLAGS)

