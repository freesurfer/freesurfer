## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_matrix_multiply
mri_matrix_multiply_SOURCES=mri_matrix_multiply.c
mri_matrix_multiply_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_matrix_multiply_LDFLAGS=$(OS_LDFLAGS)

