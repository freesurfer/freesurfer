##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom

bin_PROGRAMS = mri_binarize
mri_binarize_SOURCES=mri_binarize.c
mri_binarize_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_binarize_LDFLAGS=$(OS_LDFLAGS)

EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
