## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_LDFLAGS=

bin_PROGRAMS = mris_convert
mris_convert_SOURCES=mris_convert.c
mris_convert_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mris_convert_LDFLAGS=$(OS_LDFLAGS)

# trick to get test data into the build directory (mris_expand testdat is used)
foo_SOURCES=
foo:
	if ! test -f "$(top_builddir)/mris_expand/testdata.tar.gz"; then \
		cp -v $(top_srcdir)/mris_expand/testdata.tar.gz \
					$(top_builddir)/mris_expand/; fi;

check_PROGRAMS=foo

TESTS=test_mris_convert

EXTRA_DIST=test_mris_convert

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
