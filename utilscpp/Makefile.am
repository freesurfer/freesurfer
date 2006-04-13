#
# Makefile.am for libutilscpp.a - specify in addition to the standard
#
AM_CPPFLAGS=$(X_CFLAGS) \
	-I$(top_srcdir)/include \
	-I$(top_srcdir)/include/dicom

noinst_LIBRARIES= libutilscpp.a
libutilscpp_a_SOURCES = mri_fastmarching.cpp mris_fastmarching.cpp globals.cpp face.cpp vertex.cpp loop.cpp fastloop.cpp surface.cpp patchdisk.cpp 
 
SUBDIRS=test


# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
