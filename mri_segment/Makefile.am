## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_LDFLAGS=

bin_PROGRAMS = mri_segment

mri_segment_SOURCES=mri_segment.c mri_segment.help.xml.h
mri_segment_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_segment_LDFLAGS=$(OS_LDFLAGS)

foodir=$(prefix)/docs/xml
foo_DATA=mri_segment.help.xml
foo2dir=$(prefix)/docs/html
foo2_DATA=mri_segment.help.xml.html

TESTS=$(top_builddir)/scripts/help_xml_validate

clean-local:
	rm -f $(BUILT_SOURCES) 

EXTRA_DIST=$(foo_DATA) $(BUILT_SOURCES)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra

BUILT_SOURCES=mri_segment.help.xml.h mri_segment.help.xml.html

