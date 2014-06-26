## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=

bin_PROGRAMS = mri_vol2surf
mri_vol2surf_SOURCES=mri_vol2surf.c
mri_vol2surf_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_vol2surf_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
