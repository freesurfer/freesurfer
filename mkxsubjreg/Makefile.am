## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mkxsubjreg
mkxsubjreg_SOURCES=mkxsubjreg.c
mkxsubjreg_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mkxsubjreg_LDFLAGS=$(OS_LDFLAGS)

