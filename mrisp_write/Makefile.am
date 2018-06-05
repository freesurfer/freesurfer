## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_LDFLAGS=

bin_PROGRAMS = mrisp_write

mrisp_write_SOURCES=mrisp_write.c mrisp_write.help.xml.h
mrisp_write_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mrisp_write_LDFLAGS=$(OS_LDFLAGS)

foodir=$(prefix)/docs/xml
foo_DATA=mrisp_write.help.xml
foo2dir=$(prefix)/docs/html
foo2_DATA=mrisp_write.help.xml.html

TESTS=$(top_builddir)/scripts/help_xml_validate

clean-local:
	rm -f $(BUILT_SOURCES) 

EXTRA_DIST=$(foo_DATA) $(BUILT_SOURCES)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra

BUILT_SOURCES=mrisp_write.help.xml.h mrisp_write.help.xml.html

