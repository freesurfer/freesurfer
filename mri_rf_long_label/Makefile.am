##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = mri_rf_long_label
#mri_rf_long_label_SOURCES=mri_rf_long_label.c mri_rf_long_label.help.xml.h
mri_rf_long_label_SOURCES=mri_rf_long_label.c 
mri_rf_long_label_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_rf_long_label_LDFLAGS=$(OS_LDFLAGS)

#foodir=$(prefix)/docs/xml
#foo_DATA=mri_rf_long_label.help.xml

#TESTS=$(top_builddir)/scripts/help_xml_validate
TESTS=

#EXTRA_DIST=$(foo_DATA) $(BUILT_SOURCES)
EXTRA_DIST=$(BUILT_SOURCES)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra

#BUILT_SOURCES=mri_rf_long_label.help.xml.h
