## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_probe_ima
mri_probe_ima_SOURCES=mri_probe_ima.c
mri_probe_ima_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_probe_ima_LDFLAGS=$(OS_LDFLAGS)

