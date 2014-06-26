## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=

bin_PROGRAMS = mri_xvolavg
mri_xvolavg_SOURCES=mri_xvolavg.c
mri_xvolavg_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_xvolavg_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
