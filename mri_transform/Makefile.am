## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_transform
mri_transform_SOURCES=mri_transform.c
mri_transform_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_transform_LDFLAGS=$(OS_LDFLAGS)

