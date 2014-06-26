## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=

bin_PROGRAMS = mri_strip_subject_info
mri_strip_subject_info_SOURCES=mri_strip_subject_info.c
mri_strip_subject_info_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_strip_subject_info_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
