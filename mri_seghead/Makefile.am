## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_seghead
mri_seghead_SOURCES=mri_seghead.c
mri_seghead_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_seghead_LDFLAGS=$(OS_LDFLAGS)

