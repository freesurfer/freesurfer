## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mris_transform
mris_transform_SOURCES=mris_transform.c
mris_transform_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_transform_LDFLAGS=$(OS_LDFLAGS)

