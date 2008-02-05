#
# Makefile.am for libutils.a
#

AM_CPPFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=

noinst_LIBRARIES= librgb.a
librgb_a_SOURCES = close.c   name.c  rdwr.c  row.c \
	filbuf.c  open.c  rgbutils.c flsbuf.c  pix.c   rle.c

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
