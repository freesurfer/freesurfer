#
# Makefile for fsgdf
#
# Makefile.am specify i addtion to the standard
AM_CPPFLAGS=-I$(top_srcdir)/include
AM_LDFLAGS=

noinst_LIBRARIES= libfsgdf.a
libfsgdf_a_SOURCES = fsgdf.c

check_PROGRAMS=test
test_SOURCES=test.c
test_LDADD=libfsgdf.a $(addprefix $(top_srcdir)/, $(LIBS_MGH))
TESTS=test
