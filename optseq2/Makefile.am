## 
## Makefile.am 
##
AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom $(GSL_CFLAGS)
AM_LDFLAGS=

bin_PROGRAMS = optseq2
optseq2_SOURCES=optseq2.c
optseq2_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
optseq2_LDFLAGS=$(OS_LDFLAGS)

TESTS=test_optseq2.csh

EXTRA_DIST=test_optseq2.csh test_data.tar.gz

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
