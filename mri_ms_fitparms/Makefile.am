## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_ms_fitparms
mri_ms_fitparms_SOURCES=mri_ms_fitparms.c
mri_ms_fitparms_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_ms_fitparms_LDFLAGS=$(OS_LDFLAGS)

