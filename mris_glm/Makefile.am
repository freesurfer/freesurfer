## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mris_glm
mris_glm_SOURCES=mris_glm.c
mris_glm_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH)) ../fsgdf/libtclfsgdf$(DYNAMIC_LIB_EXT)
mris_glm_LDFLAGS=

