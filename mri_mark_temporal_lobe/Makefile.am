## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_mark_temporal_lobe
mri_mark_temporal_lobe_SOURCES=mri_mark_temporal_lobe.c
mri_mark_temporal_lobe_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_mark_temporal_lobe_LDFLAGS=$(OS_LDFLAGS)

