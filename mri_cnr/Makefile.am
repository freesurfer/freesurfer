## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=

bin_PROGRAMS = mri_cnr
mri_cnr_SOURCES=mri_cnr.c
mri_cnr_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_cnr_LDFLAGS=$(OS_LDFLAGS)

