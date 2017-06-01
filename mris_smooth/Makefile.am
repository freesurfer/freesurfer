## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_LDFLAGS=

bin_PROGRAMS = mris_smooth

mris_smooth_SOURCES=mris_smooth.c
mris_smooth_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mris_smooth_LDFLAGS=$(OS_LDFLAGS)

foodir=$(prefix)/docs/xml
foo_DATA=mris_smooth.help.xml
foo2dir=$(prefix)/docs/html
foo2_DATA=mris_smooth.help.xml.html

# trick to get test data into the build directory
foo_SOURCES=
foo:
	if ! test -f "$(top_builddir)/mris_registetr/testdata.tar.gz"; then \
		cp -v $(top_srcdir)/mris_registetr/testdata.tar.gz \
					$(top_builddir)/mris_registetr/; fi;

TESTS=$(top_builddir)/scripts/help_xml_validate test_fs_posix

EXTRA_DIST=$(foo_DATA) $(BUILT_SOURCES)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra

BUILT_SOURCES=mris_smooth.help.xml.h

clean-local:
	-rm *.help.xml.h *.help.xml.html
