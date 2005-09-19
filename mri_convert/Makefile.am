## 
## Makefile.am 
##
AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom

bin_PROGRAMS = mri_convert
mri_convert_SOURCES=mri_convert.c
mri_convert_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_convert_LDFLAGS=$(OS_LDFLAGS)

