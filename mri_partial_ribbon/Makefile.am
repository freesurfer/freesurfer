## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=

bin_PROGRAMS = mri_partial_ribbon
mri_partial_ribbon_SOURCES=mri_partial_ribbon.c
mri_partial_ribbon_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_partial_ribbon_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
