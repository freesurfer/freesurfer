## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mris_anatomical_stats
mris_anatomical_stats_SOURCES=mris_anatomical_stats.c
mris_anatomical_stats_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_anatomical_stats_LDFLAGS=$(OS_LDFLAGS)

