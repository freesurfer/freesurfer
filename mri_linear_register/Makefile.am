## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_linear_register
mri_linear_register_SOURCES=mri_linear_register.c
mri_linear_register_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_linear_register_LDFLAGS=$(OS_LDFLAGS)

