## 
# Makefile.am
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_LDFLAGS=

bin_PROGRAMS = mri_add_xform_to_header
mri_add_xform_to_header_SOURCES=mri_add_xform_to_header.c
mri_add_xform_to_header_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_add_xform_to_header_LDFLAGS=$(OS_LDFLAGS)

foodir=$(prefix)/docs/xml
foo_DATA=mri_add_xform_to_header.help.xml

TESTS=$(top_builddir)/scripts/help_xml_validate

EXTRA_DIST=$(foo_DATA)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra

BUILT_SOURCES=mri_add_xform_to_header.help.xml.h

