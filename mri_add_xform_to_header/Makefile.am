## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_add_xform_to_header
mri_add_xform_to_header_SOURCES=mri_add_xform_to_header.c
mri_add_xform_to_header_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_add_xform_to_header_LDFLAGS=$(OS_LDFLAGS)

