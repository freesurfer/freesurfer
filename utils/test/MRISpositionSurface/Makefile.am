## 
## Makefile.am 
##

AM_CPPFLAGS=-I$(top_srcdir)/include
AM_LDFLAGS=

check_PROGRAMS = test_MRISpositionSurface

TESTS=run_test_MRISpositionSurface

test_MRISpositionSurface_SOURCES=test_MRISpositionSurface.cpp
test_MRISpositionSurface_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
test_MRISpositionSurface_LDFLAGS= $(OS_LDFLAGS)

EXTRA_DIST=run_test_MRISpositionSurface testdata.tar.gz

EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra

clean-local:
	rm -f *.o
