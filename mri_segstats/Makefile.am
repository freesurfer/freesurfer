## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom

bin_PROGRAMS = mri_segstats
mri_segstats_SOURCES=mri_segstats.c
mri_segstats_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_segstats_LDFLAGS=$(OS_LDFLAGS)

# trick to get test data into the build directory
foo_SOURCES=
foo:
	if ! test -f "$(top_builddir)/mri_segstats/testdata.tar.gz"; then \
		cp -v $(top_srcdir)/mri_segstats/testdata.tar.gz \
					$(top_builddir)/mri_segstats/; fi;

check_PROGRAMS=foo

TESTS=test_mri_segstats

EXTRA_DIST=test_mri_segstats testdata.tar.gz

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
