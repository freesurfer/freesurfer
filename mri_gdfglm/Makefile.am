## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_gdfglm
mri_gdfglm_SOURCES=mri_gdfglm.c
mri_gdfglm_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH)) $(top_srcdir)/fsgdf/libtclfsgdf$(DYNAMIC_LIB_EXT)
mri_gdfglm_LDFLAGS=$(OS_LDFLAGS)

