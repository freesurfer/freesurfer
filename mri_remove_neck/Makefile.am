## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_remove_neck
mri_remove_neck_SOURCES=mri_remove_neck.c
mri_remove_neck_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_remove_neck_LDFLAGS=$(OS_LDFLAGS)

