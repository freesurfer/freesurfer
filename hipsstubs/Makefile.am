#
# Makefile.am
#

AM_CPPFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=

noinst_LIBRARIES= libhipsstubs.a
libhipsstubs_a_SOURCES = hipsrepl.c  hipsstubs.c

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
