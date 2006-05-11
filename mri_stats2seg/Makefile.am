##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom

bin_PROGRAMS = mri_stats2seg
mri_stats2seg_SOURCES=mri_stats2seg.c
mri_stats2seg_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_stats2seg_LDFLAGS=$(OS_LDFLAGS)

EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
