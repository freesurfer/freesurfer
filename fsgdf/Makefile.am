#
# Makefile for fsgdf - specify in addition to the standard
#
AM_CPPFLAGS=-I$(top_srcdir)/include
AM_LDFLAGS=

noinst_LIBRARIES= libfsgdf.a
libfsgdf_a_SOURCES = fsgdf.c

check_PROGRAMS=test
test_SOURCES=test.c fsgdf_wrap.c
test_LDFLAGS=$(LDFLAGS_TCL)
test_LDADD=libfsgdf.a $(addprefix $(top_srcdir)/, $(LIBS_MGH)) $(LIBS_TCL)
TESTS=test


# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
