## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_make_register
mri_make_register_SOURCES=mri_make_register.c
mri_make_register_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_make_register_LDFLAGS=$(OS_LDFLAGS)

