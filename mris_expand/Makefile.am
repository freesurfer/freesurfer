## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = mris_expand
mris_expand_SOURCES=mris_expand.c
mris_expand_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_expand_LDFLAGS=$(OS_LDFLAGS)

TESTS=test_mris_expand

EXTRA_DIST=testdata.tar.gz test_mris_expand

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
