##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS=mri_warp_concat
mri_warp_concat_SOURCES=mri_warp_concat.cpp mri_warp_concat.help.xml.h
mri_warp_concat_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_warp_concat_LDFLAGS=$(OS_LDFLAGS)

foodir=$(prefix)/docs/xml
foo_DATA=mri_warp_concat.help.xml

TESTS=$(top_builddir)/scripts/help_xml_validate

EXTRA_DIST=$(foo_DATA) $(BUILT_SOURCES)

EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra

BUILT_SOURCES=mri_warp_concat.help.xml.h

clean-local:
	rm -f $(BUILT_SOURCES)
