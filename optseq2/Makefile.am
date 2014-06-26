## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_LDFLAGS=

bin_PROGRAMS = optseq2
optseq2_SOURCES=optseq2.c
optseq2_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
optseq2_LDFLAGS=$(OS_LDFLAGS)

# trick to get test data into the build directory
foo_SOURCES=
foo:
	if ! test -f "$(top_builddir)/optseq2/test_data.tar.gz"; then \
		cp -v $(top_srcdir)/optseq2/test_data.tar.gz \
					$(top_builddir)/optseq2/; fi;

check_PROGRAMS=foo

TESTS=test_optseq2.csh

EXTRA_DIST=test_optseq2.csh test_data.tar.gz

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
