## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=

bin_PROGRAMS = mris_find_flat_regions
mris_find_flat_regions_SOURCES=mris_find_flat_regions.c
mris_find_flat_regions_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_find_flat_regions_LDFLAGS=$(OS_LDFLAGS)



# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
