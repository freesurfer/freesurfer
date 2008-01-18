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
AM_CPPFLAGS=-I$(top_srcdir)/include
endif
AM_LDFLAGS=

noinst_LIBRARIES= libfsgdf.a
libfsgdf_a_SOURCES = fsgdf.c

if ENABLE_TCLTK_APPS
check_PROGRAMS=test
test_SOURCES=test.c fsgdf_wrap.c
test_LDFLAGS=$(LDFLAGS_TCL)
test_LDADD=libfsgdf.a $(addprefix $(top_srcdir)/, $(LIBS_MGH)) $(LIBS_TCL)
TESTS=setup_test_data test
endif

EXTRA_DIST=setup_test_data test_data.tar.gz fsgdf.i test.tcl

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
