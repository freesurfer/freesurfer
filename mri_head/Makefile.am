## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_head
mri_head_SOURCES=mri_head.c
mri_head_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_head_LDFLAGS=$(OS_LDFLAGS)

