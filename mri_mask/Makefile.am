## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_mask
mri_mask_SOURCES=mri_mask.c
mri_mask_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_mask_LDFLAGS=$(OS_LDFLAGS)

