## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mris_morph_stats
mris_morph_stats_SOURCES=mris_morph_stats.c
mris_morph_stats_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_morph_stats_LDFLAGS=$(OS_LDFLAGS)

