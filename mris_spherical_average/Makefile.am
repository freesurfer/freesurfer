## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mris_spherical_average
mris_spherical_average_SOURCES=mris_spherical_average.c
mris_spherical_average_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_spherical_average_LDFLAGS=$(OS_LDFLAGS)

