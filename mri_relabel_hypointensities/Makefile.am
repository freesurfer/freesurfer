## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=

bin_PROGRAMS = mri_relabel_hypointensities
mri_relabel_hypointensities_SOURCES=mri_relabel_hypointensities.c
mri_relabel_hypointensities_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_relabel_hypointensities_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
