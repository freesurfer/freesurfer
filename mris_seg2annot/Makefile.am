##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom

bin_PROGRAMS = mris_seg2annot
mris_seg2annot_SOURCES=mris_seg2annot.c
mris_seg2annot_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mris_seg2annot_LDFLAGS=$(OS_LDFLAGS)

EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
