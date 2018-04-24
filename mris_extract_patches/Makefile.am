## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = mris_extract_patches
mris_extract_patches_SOURCES=mris_extract_patches.c
mris_extract_patches_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mris_extract_patches_LDFLAGS=$(OS_LDFLAGS)

# trick to get test data into the build directory
foo_SOURCES=
foo:
	if ! test -f "$(top_builddir)/mris_extract_patches/testdata.tar.gz"; then \
		cp -v $(top_srcdir)/mris_extract_patches/testdata.tar.gz \
					$(top_builddir)/mris_extract_patches/; fi;

check_PROGRAMS=foo

TESTS=test_mris_extract_patches

EXTRA_DIST=testdata.tar.gz test_mris_extract_patches

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
