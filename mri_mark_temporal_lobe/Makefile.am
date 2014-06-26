## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=

bin_PROGRAMS = mri_mark_temporal_lobe
mri_mark_temporal_lobe_SOURCES=mri_mark_temporal_lobe.c
mri_mark_temporal_lobe_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_mark_temporal_lobe_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
