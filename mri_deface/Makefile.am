## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=

bin_PROGRAMS = mri_deface
mri_deface_SOURCES=mri_deface.c
mri_deface_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_deface_LDFLAGS=$(OS_LDFLAGS)

