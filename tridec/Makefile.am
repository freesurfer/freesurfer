## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = tridec
tridec_SOURCES=tridec.c
tridec_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
tridec_LDFLAGS=$(OS_LDFLAGS)

