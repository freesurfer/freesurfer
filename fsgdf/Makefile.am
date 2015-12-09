#
# Makefile for fsgdf
#

#
# Note: fsgdf_wrap.c is created using swig:
#   swig -tcl fsgdf.i
# and then all the unused routines, as indicated by the compiler errors,
# are removed.  the callocs were also increased to size 2048 for safety.
# fsgdf_wrap.c is built by tksurfer and qdec.
#

if ENABLE_TCLTK_APPS
AM_CPPFLAGS=-I$(top_srcdir)/include $(TCL_CFLAGS)
else
AM_CPPFLAGS=-I$(top_srcdir)/include $(TCL_CFLAGS)
endif
AM_LDFLAGS=

noinst_LIBRARIES= libfsgdf.a
if ENABLE_TCLTK_APPS
libfsgdf_a_SOURCES = fsgdf.c fsgdf_wrap.c FsgdfPlot.cxx
else
libfsgdf_a_SOURCES = fsgdf.c fsgdf_wrap.c
endif

if ENABLE_TCLTK_APPS
# trick to get test data into the build directory
foo_SOURCES=
foo:
	if ! test -f "$(top_builddir)/fsgdf/test_data.tar.gz"; then \
		cp -v $(top_srcdir)/fsgdf/test_data.tar.gz \
					$(top_builddir)/fsgdf/; fi;
	if ! test -f "$(top_builddir)/fsgdf/test.tcl"; then \
		cp -v $(top_srcdir)/fsgdf/test.tcl \
					$(top_builddir)/fsgdf/; fi;

check_PROGRAMS=foo test
test_SOURCES=test.c fsgdf_wrap.c
test_LDFLAGS=$(LDFLAGS_TCL)
test_LDADD=libfsgdf.a $(addprefix $(top_builddir)/, $(LIBS_MGH)) $(LIBS_TCL)
TESTS=setup_test_data test
endif

EXTRA_DIST=setup_test_data test_data.tar.gz fsgdf.i test.tcl FsgdfPlot.cxx

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
