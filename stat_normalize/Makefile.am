## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = stat_normalize
stat_normalize_SOURCES=stat_normalize.c
stat_normalize_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
stat_normalize_LDFLAGS=$(OS_LDFLAGS)

