## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_paint
mri_paint_SOURCES=mri_paint.c
mri_paint_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_paint_LDFLAGS=$(OS_LDFLAGS)

