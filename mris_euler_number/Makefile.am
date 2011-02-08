## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_LDFLAGS=

bin_PROGRAMS = mris_euler_number

mris_euler_number_SOURCES=mris_euler_number.c mris_euler_number.help.xml.h
mris_euler_number_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mris_euler_number_LDFLAGS=$(OS_LDFLAGS)

foodir=$(prefix)/docs/xml
foo_DATA=mris_euler_number.help.xml

TESTS=$(top_builddir)/scripts/help_xml_validate

EXTRA_DIST=$(foo_DATA) $(BUILT_SOURCES)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra

BUILT_SOURCES=mris_euler_number.help.xml.h

