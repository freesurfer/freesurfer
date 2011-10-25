##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom

bin_PROGRAMS = mri_extract_largest_CC 
mri_extract_largest_CC_SOURCES=mri_extract_largest_CC.c 
mri_extract_largest_CC_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_extract_largest_CC_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
