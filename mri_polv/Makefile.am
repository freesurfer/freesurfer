## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_polv
mri_polv_SOURCES=mri_polv.c
mri_polv_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_polv_LDFLAGS=$(OS_LDFLAGS)

