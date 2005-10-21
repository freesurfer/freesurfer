## 
# Makefile.am - specify in addition to the standard
##

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=

bin_PROGRAMS = mri_add_xform_to_header
mri_add_xform_to_header_SOURCES=mri_add_xform_to_header.c
mri_add_xform_to_header_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_add_xform_to_header_LDFLAGS=$(OS_LDFLAGS)



# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
