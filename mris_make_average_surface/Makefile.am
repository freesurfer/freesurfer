## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mris_make_average_surface
mris_make_average_surface_SOURCES=mris_make_average_surface.c
mris_make_average_surface_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_make_average_surface_LDFLAGS=$(OS_LDFLAGS)

