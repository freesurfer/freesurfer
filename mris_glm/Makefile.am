## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mris_glm
mris_glm_SOURCES=mris_glm.c
mris_glm_LDADD= $(top_srcdir)/fsgdf/libfsgdf.a $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_glm_LDFLAGS=

