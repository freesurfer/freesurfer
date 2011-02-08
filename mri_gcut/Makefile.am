AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = mri_gcut

mri_gcut_SOURCES=mri_gcut.cpp mri_gcut.help.xml.h
mri_gcut_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_gcut_LDFLAGS=$(OS_LDFLAGS)

foodir=$(prefix)/docs/xml
foo_DATA=mri_gcut.help.xml

TESTS=$(top_builddir)/scripts/help_xml_validate

EXTRA_DIST = gcut.h graphcut.h graphcut.cpp pre_pro.cpp \
	$(foo_DATA) $(BUILT_SOURCES)

EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra

BUILT_SOURCES=mri_gcut.help.xml.h

