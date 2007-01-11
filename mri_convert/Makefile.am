## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom

bin_PROGRAMS = mri_convert
mri_convert_SOURCES=mri_convert.c
mri_convert_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_convert_LDFLAGS=$(OS_LDFLAGS)

TESTS=test_mri_convert

EXTRA_DIST=test_mri_convert testdata.tar.gz

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
