## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_LDFLAGS=

bin_PROGRAMS = mris_ca_label

mris_ca_label_SOURCES=mris_ca_label.c mris_ca_label.help.xml.h
mris_ca_label_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mris_ca_label_LDFLAGS=$(OS_LDFLAGS)

foodir=$(prefix)/docs/xml
foo_DATA=mris_ca_label.help.xml

# trick to get test data into the build directory
foo_SOURCES=
foo:
	if ! test -f "$(top_builddir)/mris_ca_label/testdata.tar.gz"; then \
		cp -v $(top_srcdir)/mris_ca_label/testdata.tar.gz \
					$(top_builddir)/mris_ca_label/; fi;
	if ! test -f "$(top_builddir)/distribution/average/lh.curvature.buckner40.filled.desikan_killiany.2010-03-25.gcs"; then \
		cp -v $(top_srcdir)/distribution/average/*.gcs \
					$(top_builddir)/distribution/average/; fi;

check_PROGRAMS=foo

TESTS=test_mris_ca_label $(top_builddir)/scripts/help_xml_validate

EXTRA_DIST=test_mris_ca_label testdata.tar.gz $(foo_DATA) $(BUILT_SOURCES)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra

BUILT_SOURCES=mris_ca_label.help.xml.h

