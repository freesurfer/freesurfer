## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_edit_segmentation_with_surfaces
mri_edit_segmentation_with_surfaces_SOURCES=mri_edit_segmentation_with_surfaces.c
mri_edit_segmentation_with_surfaces_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_edit_segmentation_with_surfaces_LDFLAGS=$(OS_LDFLAGS)

