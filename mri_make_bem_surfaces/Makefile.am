## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_LDFLAGS=

bin_PROGRAMS = mri_make_bem_surfaces
mri_make_bem_surfaces_SOURCES=mri_make_bem_surfaces.c ../utils/machine.c
# cyclic dependency between rgb and utils on swapInt etc.
mri_make_bem_surfaces_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_make_bem_surfaces_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
