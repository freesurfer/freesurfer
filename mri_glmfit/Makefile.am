## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = mri_glmfit
mri_glmfit_SOURCES=mri_glmfit.c $(top_srcdir)/fsgdf/fsgdf.c
mri_glmfit_LDADD=$(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_glmfit_LDFLAGS=$(OS_LDFLAGS)

# trick to get test data into the build directory
foo_SOURCES=
foo:
	if ! test -f "$(top_builddir)/mri_glmfit/teststats.tar.gz"; then \
		cp -v $(top_srcdir)/mri_glmfit/teststats.tar.gz \
					$(top_builddir)/mri_glmfit/; fi;

check_PROGRAMS=foo

TESTS=test_mri_glmfit

EXTRA_DIST=test_mri_glmfit teststats.tar.gz

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
