## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mris_annot_to_segmentation
mris_annot_to_segmentation_SOURCES=mris_annot_to_segmentation.c
mris_annot_to_segmentation_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_annot_to_segmentation_LDFLAGS=$(OS_LDFLAGS)

