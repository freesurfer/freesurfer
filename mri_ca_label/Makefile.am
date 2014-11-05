## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_LDFLAGS=

bin_PROGRAMS = mri_ca_label

mri_ca_label_SOURCES=mri_ca_label.c mri_ca_label.help.xml.h
mri_ca_label_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_ca_label_LDFLAGS=$(OS_LDFLAGS)

foodir=$(prefix)/docs/xml
foo_DATA=mri_ca_label.help.xml
foo2dir=$(prefix)/docs/html
foo2_DATA=mri_ca_label.help.xml.html

# this one weird trick to get test data into the build directory
# that doctors don't want you to know about
foo_SOURCES=
foo:
	if ! test -f "$(top_builddir)/mri_ca_label/testdata.tar.gz"; then \
		cp -v $(top_srcdir)/mri_ca_label/testdata.tar.gz \
					$(top_builddir)/mri_ca_label/; fi;

TESTS=test_mri_ca_label $(top_builddir)/scripts/help_xml_validate

clean-local:
	rm -f $(BUILT_SOURCES)

EXTRA_DIST=$(foo_DATA) $(BUILT_SOURCES) test_mri_ca_label

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra

BUILT_SOURCES=mri_ca_label.help.xml.h mri_ca_label.help.xml.html
