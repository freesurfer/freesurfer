## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mris_translate_annotation
mris_translate_annotation_SOURCES=mris_translate_annotation.c
mris_translate_annotation_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_translate_annotation_LDFLAGS=$(OS_LDFLAGS)

