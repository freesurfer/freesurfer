## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mris_surface_to_vol_distances
mris_surface_to_vol_distances_SOURCES=mris_surface_to_vol_distances.c
mris_surface_to_vol_distances_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_surface_to_vol_distances_LDFLAGS=$(OS_LDFLAGS)

