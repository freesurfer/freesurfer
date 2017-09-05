##
## Makefile.am
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

BUILT_SOURCES=mris_watershed.help.xml.h mris_watershed.help.xml.html

foodir=$(prefix)/docs/xml
foo_DATA=mris_watershed.help.xml
foo2dir=$(prefix)/docs/html
foo2_DATA=mris_watershed.help.xml.html

bin_PROGRAMS = mris_watershed
mris_watershed_SOURCES=mris_watershed.c
mris_watershed_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mris_watershed_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra

clean-local:
	rm -f $(BUILT_SOURCES)
