## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mrisp_paint
mrisp_paint_SOURCES=mrisp_paint.c
mrisp_paint_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mrisp_paint_LDFLAGS=$(OS_LDFLAGS)

