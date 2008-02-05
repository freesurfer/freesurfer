#
# Makefile.am
#

AM_CPPFLAGS=-I$(top_srcdir)/include
AM_LDFLAGS=

noinst_LIBRARIES= liblog.a
liblog_a_SOURCES = congraph.c  h_logz.c   map.c \
	connect.c   lmedian.c  runfuncs.c

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
