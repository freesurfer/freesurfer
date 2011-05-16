## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_LDFLAGS=

bin_PROGRAMS = mri_ms_fitparms
mri_ms_fitparms_SOURCES=mri_ms_fitparms.c mri_ms_fitparms.help.xml.h
mri_ms_fitparms_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_ms_fitparms_LDFLAGS=$(OS_LDFLAGS)

foodir=$(prefix)/docs/xml
foo_DATA=mri_ms_fitparms.help.xml

TESTS=$(top_builddir)/scripts/help_xml_validate

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra

BUILT_SOURCES=mri_ms_fitparms.help.xml.h