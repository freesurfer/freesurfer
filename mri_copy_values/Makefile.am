## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_copy_values
mri_copy_values_SOURCES=mri_copy_values.c
mri_copy_values_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_copy_values_LDFLAGS=$(OS_LDFLAGS)

