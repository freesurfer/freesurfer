## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_vol2roi
mri_vol2roi_SOURCES=mri_vol2roi.c
mri_vol2roi_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_vol2roi_LDFLAGS=$(OS_LDFLAGS)

