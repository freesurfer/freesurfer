#
# Makefile.am for libutilscpp.a
#

AM_CPPFLAGS=-std=c++0x $(X_CFLAGS) \
	-I$(top_srcdir)/include \
	-I$(top_srcdir)/include/dicom

SUBDIRS=test

noinst_LIBRARIES= libutilscpp.a

libutilscpp_a_SOURCES = argparse.cpp log.cpp mris_fastmarching.cpp globals.cpp face.cpp \
	vertex.cpp loop.cpp fastloop.cpp surface.cpp patchdisk.cpp \
	mris_topology.cpp segment.cpp patchdisk.cpp

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
