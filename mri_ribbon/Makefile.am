## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_ribbon
mri_ribbon_SOURCES=mri_ribbon.c
mri_ribbon_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_ribbon_LDFLAGS=$(OS_LDFLAGS)

