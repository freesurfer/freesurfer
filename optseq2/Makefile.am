## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = optseq2
optseq2_SOURCES=optseq2.c
optseq2_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
optseq2_LDFLAGS=$(OS_LDFLAGS)
