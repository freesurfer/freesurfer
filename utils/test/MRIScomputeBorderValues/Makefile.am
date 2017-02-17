## 
## Makefile.am 
##

AM_CPPFLAGS=-I$(top_srcdir)/include
AM_LDFLAGS=

check_PROGRAMS = test_MRIScomputeBorderValues

TESTS=run_test_MRIScomputeBorderValues

test_MRIScomputeBorderValues_SOURCES=test_MRIScomputeBorderValues.cpp
test_MRIScomputeBorderValues_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
test_MRIScomputeBorderValues_LDFLAGS= $(OS_LDFLAGS)

EXTRA_DIST=run_test_MRIScomputeBorderValues testdata.tar.gz

EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra

clean-local:
	rm -f *.o
