## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_strip_nonwhite
mri_strip_nonwhite_SOURCES=mri_strip_nonwhite.c
mri_strip_nonwhite_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_strip_nonwhite_LDFLAGS=$(OS_LDFLAGS)

