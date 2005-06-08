#
# Makefile.am for libutilscpp.a - specify in addition to the standard
#
AM_CPPFLAGS=$(X_CFLAGS) \
	-I$(top_srcdir)/include \
	-I$(top_srcdir)/include/dicom

noinst_LIBRARIES= libutilscpp.a
libutilscpp_a_SOURCES = mri_fastmarching.cpp mris_fastmarching.cpp

SUBDIRS=test
