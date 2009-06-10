##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = mri_concatenate_lta
mri_concatenate_lta_SOURCES=mri_concatenate_lta.c
mri_concatenate_lta_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_concatenate_lta_LDFLAGS=$(OS_LDFLAGS)

# trick to get test data into the build directory during make distcheck
foo_SOURCES=
foo:
	if ! test -f "$(top_builddir)/mri_concatenate_lta/testdata.tar.gz"; then \
		cp -vf $(top_srcdir)/mri_concatenate_lta/testdata.tar.gz \
					$(top_builddir)/mri_concatenate_lta/; fi;

check_PROGRAMS=foo

TESTS=test_mri_concatenate_lta

EXTRA_DIST=test_mri_concatenate_lta testdata.tar.gz

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
