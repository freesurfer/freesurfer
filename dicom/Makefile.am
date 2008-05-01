#
# Makefile.am
#

AM_CPPFLAGS=-I$(top_srcdir)/include \
	-I$(top_srcdir)/include/dicom

AM_LDFLAGS=

noinst_LIBRARIES= libdicom.a

libdicom_a_SOURCES = condition.c dcmcond.c lst.c \
	ctnthread.c  dcmdict.c dcm.c dcmsupport.c

EXTRA_DIST=

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
