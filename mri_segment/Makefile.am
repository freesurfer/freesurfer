## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_segment
mri_segment_SOURCES=mri_segment.c
mri_segment_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_segment_LDFLAGS=$(OS_LDFLAGS)

