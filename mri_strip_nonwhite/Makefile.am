## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=

bin_PROGRAMS = mri_strip_nonwhite
mri_strip_nonwhite_SOURCES=mri_strip_nonwhite.c
mri_strip_nonwhite_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_strip_nonwhite_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
