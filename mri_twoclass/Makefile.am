## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_twoclass
mri_twoclass_SOURCES=mri_twoclass.c
mri_twoclass_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_twoclass_LDFLAGS=$(OS_LDFLAGS)

