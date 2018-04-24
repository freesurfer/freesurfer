## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=

bin_PROGRAMS = mris_transform
mris_transform_SOURCES=mris_transform.c
mris_transform_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mris_transform_LDFLAGS=$(OS_LDFLAGS)

TESTS=test

EXTRA_DIST=test testdata.tar.gz $(foo_DATA) $(BUILT_SOURCES)

EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
