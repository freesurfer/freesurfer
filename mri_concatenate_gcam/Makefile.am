##
## Makefile.am
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS=mri_concatenate_gcam
mri_concatenate_gcam_SOURCES=mri_concatenate_gcam.cpp mri_concatenate_gcam.help.xml.h
mri_concatenate_gcam_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_concatenate_gcam_LDFLAGS=$(OS_LDFLAGS)

foodir=$(prefix)/docs/xml
foo_DATA=mri_concatenate_gcam.help.xml

TESTS=$(top_builddir)/scripts/help_xml_validate

EXTRA_DIST=$(foo_DATA) $(BUILT_SOURCES)

EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra

BUILT_SOURCES=mri_concatenate_gcam.help.xml.h

clean-local:
rm -f $(BUILT_SOURCES)
