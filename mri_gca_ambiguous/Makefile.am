## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_gca_ambiguous
mri_gca_ambiguous_SOURCES=mri_gca_ambiguous.c
mri_gca_ambiguous_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_gca_ambiguous_LDFLAGS=$(OS_LDFLAGS)

