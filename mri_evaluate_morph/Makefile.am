## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_evaluate_morph
mri_evaluate_morph_SOURCES=mri_evaluate_morph.c
mri_evaluate_morph_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_evaluate_morph_LDFLAGS=$(OS_LDFLAGS)

