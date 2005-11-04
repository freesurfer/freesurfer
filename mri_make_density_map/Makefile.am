## 
## Makefile.am 
##
AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom

bin_PROGRAMS = mri_make_density_map
mri_make_density_map_SOURCES=mri_make_density_map.c
mri_make_density_map_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_make_density_map_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
