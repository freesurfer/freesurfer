## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_strip_subject_info
mri_strip_subject_info_SOURCES=mri_strip_subject_info.c
mri_strip_subject_info_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_strip_subject_info_LDFLAGS=$(OS_LDFLAGS)

