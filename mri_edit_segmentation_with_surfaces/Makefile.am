## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_LDFLAGS=

bin_PROGRAMS = mri_edit_segmentation_with_surfaces
mri_edit_segmentation_with_surfaces_SOURCES=mri_edit_segmentation_with_surfaces.c
mri_edit_segmentation_with_surfaces_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_edit_segmentation_with_surfaces_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
