## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mris_multiscale_stats
mris_multiscale_stats_SOURCES=mris_multiscale_stats.c
mris_multiscale_stats_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_multiscale_stats_LDFLAGS=$(OS_LDFLAGS)

