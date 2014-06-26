## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_LDFLAGS=

bin_PROGRAMS = mri_normalize
mri_normalize_SOURCES=mri_normalize.c mri_normalize.help.xml.h
mri_normalize_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_normalize_LDFLAGS=$(OS_LDFLAGS)

foodir=$(prefix)/docs/xml
foo_DATA=mri_normalize.help.xml
foo2dir=$(prefix)/docs/html
foo2_DATA=mri_normalize.help.xml.html

TESTS=$(top_builddir)/scripts/help_xml_validate

EXTRA_DIST=$(foo_DATA) $(BUILT_SOURCES)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra

BUILT_SOURCES=mri_normalize.help.xml.h mri_normalize.help.xml.html

clean-local:
	rm -f $(BUILT_SOURCES)
