## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_em_register
mri_em_register_SOURCES=mri_em_register.c
mri_em_register_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_em_register_LDFLAGS=$(OS_LDFLAGS)

