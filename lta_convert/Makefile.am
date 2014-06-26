##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

BUILT_SOURCES=lta_convert.help.xml.h lta_convert.help.xml.html

foodir=$(prefix)/docs/xml
foo_DATA=lta_convert.help.xml
foo2dir=$(prefix)/docs/html
foo2_DATA=lta_convert.help.xml.html

bin_PROGRAMS = lta_convert
lta_convert_SOURCES=lta_convert.cpp
lta_convert_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
lta_convert_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra

clean-local:
	rm -f $(BUILT_SOURCES)
