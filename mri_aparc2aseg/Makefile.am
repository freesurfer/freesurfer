## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = mri_aparc2aseg
mri_aparc2aseg_SOURCES=mri_aparc2aseg.c
mri_aparc2aseg_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_aparc2aseg_LDFLAGS=$(OS_LDFLAGS)

TESTS=test_mri_aparc2aseg

EXTRA_DIST=test_mri_aparc2aseg testdata.tar.gz testdata2.tar.gz

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
