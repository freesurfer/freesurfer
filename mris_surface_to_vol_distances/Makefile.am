## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=

bin_PROGRAMS = mris_surface_to_vol_distances
mris_surface_to_vol_distances_SOURCES=mris_surface_to_vol_distances.c
mris_surface_to_vol_distances_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mris_surface_to_vol_distances_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
