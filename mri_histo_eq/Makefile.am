## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_histo_eq
mri_histo_eq_SOURCES=mri_histo_eq.c
mri_histo_eq_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_histo_eq_LDFLAGS=$(OS_LDFLAGS)

