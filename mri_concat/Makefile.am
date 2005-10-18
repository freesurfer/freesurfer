## 
## Makefile.am 
##
AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom

bin_PROGRAMS = mri_concat
mri_concat_SOURCES=mri_concat.c
mri_concat_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_concat_LDFLAGS=$(OS_LDFLAGS)
