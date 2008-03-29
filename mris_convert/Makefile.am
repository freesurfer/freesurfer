## 
## Makefile.am 
##

if USE_LOCAL_EXPAT
EXPAT_THINGS=-I$(top_srcdir)/expat
else
EXPAT_THINGS=
endif

AM_CFLAGS=-I$(top_srcdir)/include $(EXPAT_THINGS)
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
	if ! test -f "$(top_builddir)/mris_convert/giftidata.tar.gz"; then \
		cp -v $(top_srcdir)/mris_convert/giftidata.tar.gz \
					$(top_builddir)/mris_convert/; fi;

check_PROGRAMS=foo

TESTS=test_mris_convert test_gifti

EXTRA_DIST=test_mris_convert test_gifti giftidata.tar.gz

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
