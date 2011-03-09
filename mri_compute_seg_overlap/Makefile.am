##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = mri_compute_seg_overlap
mri_compute_seg_overlap_SOURCES=mri_compute_seg_overlap.c \
	mri_compute_seg_overlap.help.xml.h
mri_compute_seg_overlap_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_compute_seg_overlap_LDFLAGS=$(OS_LDFLAGS)

TESTS=$(top_builddir)/scripts/help_xml_validate

foodir=$(prefix)/docs/xml
foo_DATA=mri_compute_seg_overlap.help.xml

EXTRA_DIST=$(foo_DATA)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra

BUILT_SOURCES=mri_compute_seg_overlap.help.xml.h