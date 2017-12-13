## 
## Makefile.am 
##

AM_CXXFLAGS=-std=c++11 -I$(top_srcdir)/include -D$(target_cpu)
AM_LDFLAGS=

bin_PROGRAMS=mri_segcentroids

mri_segcentroids_SOURCES=mri_segcentroids.cpp mri_segcentroids.help.xml.h
mri_segcentroids_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_segcentroids_LDFLAGS=$(OS_LDFLAGS)

foodir=$(prefix)/docs/xml
foo_DATA=mri_segcentroids.help.xml
foo2dir=$(prefix)/docs/html
foo2_DATA=mri_segcentroids.help.xml.html

TESTS=$(top_builddir)/scripts/help_xml_validate

clean-local:
	rm -f $(BUILT_SOURCES) 

EXTRA_DIST=$(foo_DATA) $(BUILT_SOURCES)

EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra

BUILT_SOURCES=mri_segcentroids.help.xml.h mri_segcentroids.help.xml.html
