## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_reduce
mri_reduce_SOURCES=mri_reduce.c
mri_reduce_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_reduce_LDFLAGS=$(OS_LDFLAGS)

