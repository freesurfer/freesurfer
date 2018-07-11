##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom

bin_PROGRAMS = mri_diff
mri_diff_SOURCES=mri_diff.c
mri_diff_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_diff_LDFLAGS=$(OS_LDFLAGS)

TESTS=test.py

EXTRA_DIST= test.py

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
