## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_xvolavg
mri_xvolavg_SOURCES=mri_xvolavg.c
mri_xvolavg_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_xvolavg_LDFLAGS=$(OS_LDFLAGS)

