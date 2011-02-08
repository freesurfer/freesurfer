## 
## Makefile.am 
##

AM_CXXFLAGS=-I$(top_srcdir)/include -D$(target_cpu)
AM_LDFLAGS=

bin_PROGRAMS = mri_watershed  mri_validate_skull_stripped

mri_watershed_SOURCES=mri_watershed.cpp mri_watershed.help.xml.h
mri_watershed_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_watershed_LDFLAGS=$(OS_LDFLAGS)

mri_validate_skull_stripped_SOURCES=mri_validate_skull_stripped.cpp
mri_validate_skull_stripped_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_validate_skull_stripped_LDFLAGS=$(OS_LDFLAGS)

foodir=$(prefix)/docs/xml
foo_DATA=mri_watershed.help.xml

SUBDIRS=brain_volume

TESTS=$(top_builddir)/scripts/help_xml_validate

EXTRA_DIST=$(foo_DATA) $(BUILT_SOURCES)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra

BUILT_SOURCES=mri_watershed.help.xml.h

