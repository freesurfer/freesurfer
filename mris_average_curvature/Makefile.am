## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mris_average_curvature
mris_average_curvature_SOURCES=mris_average_curvature.c
mris_average_curvature_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_average_curvature_LDFLAGS=$(OS_LDFLAGS)

