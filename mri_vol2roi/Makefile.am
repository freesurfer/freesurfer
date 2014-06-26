## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=

bin_PROGRAMS = mri_vol2roi
mri_vol2roi_SOURCES=mri_vol2roi.c
mri_vol2roi_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_vol2roi_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
