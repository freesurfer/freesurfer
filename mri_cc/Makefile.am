## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_cc
mri_cc_SOURCES=mri_cc.c
mri_cc_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_cc_LDFLAGS=$(OS_LDFLAGS)

