## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = mri_aparc2aseg
mri_aparc2aseg_SOURCES=mri_aparc2aseg.c
mri_aparc2aseg_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_aparc2aseg_LDFLAGS=$(OS_LDFLAGS)

foodir=$(prefix)/docs/xml
foo_DATA=mri_aparc2aseg.help.xml

# trick to get test data into the build directory
foo_SOURCES=
foo:
	if ! test -f "$(top_builddir)/mri_aparc2aseg/testdata.tar.gz"; then \
		cp -v $(top_srcdir)/mri_aparc2aseg/testdata.tar.gz \
					$(top_builddir)/mri_aparc2aseg/; fi;
	if ! test -f "$(top_builddir)/mri_aparc2aseg/testdata2.tar.gz"; then \
		cp -v $(top_srcdir)/mri_aparc2aseg/testdata2.tar.gz \
					$(top_builddir)/mri_aparc2aseg/; fi;

check_PROGRAMS=foo

TESTS=test_mri_aparc2aseg

EXTRA_DIST=test_mri_aparc2aseg testdata.tar.gz testdata2.tar.gz $(foo_DATA)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
