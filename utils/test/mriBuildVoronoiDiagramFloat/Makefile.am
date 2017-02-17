## 
## Makefile.am 
##

AM_CPPFLAGS=-I$(top_srcdir)/include
AM_LDFLAGS=

check_PROGRAMS = test_mriBuildVoronoiDiagramFloat

TESTS=test_mriBuildVoronoiDiagramFloat_withData

test_mriBuildVoronoiDiagramFloat_SOURCES=test_mriBuildVoronoiDiagramFloat.cpp
test_mriBuildVoronoiDiagramFloat_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
test_mriBuildVoronoiDiagramFloat_LDFLAGS= $(OS_LDFLAGS)

EXTRA_DIST=testdata.tar.gz

EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra

clean-local:
	rm -f *.o
