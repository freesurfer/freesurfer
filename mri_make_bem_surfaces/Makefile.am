## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_make_bem_surfaces
mri_make_bem_surfaces_SOURCES=mri_make_bem_surfaces.c ../utils/machine.c
# cyclic dependency between rgb and utils on swapInt etc.
mri_make_bem_surfaces_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_make_bem_surfaces_LDFLAGS=$(OS_LDFLAGS)

