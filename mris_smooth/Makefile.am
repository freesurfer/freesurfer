## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mris_smooth
mris_smooth_SOURCES=mris_smooth.c
mris_smooth_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_smooth_LDFLAGS=$(OS_LDFLAGS)

