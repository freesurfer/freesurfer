## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_convert_mdh
mri_convert_mdh_SOURCES=mri_convert_mdh.c
mri_convert_mdh_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_convert_mdh_LDFLAGS=$(OS_LDFLAGS)

