## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mris_rotate
mris_rotate_SOURCES=mris_rotate.c
mris_rotate_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_rotate_LDFLAGS=$(OS_LDFLAGS)

