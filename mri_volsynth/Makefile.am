## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_volsynth
mri_volsynth_SOURCES=mri_volsynth.c
mri_volsynth_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_volsynth_LDFLAGS=$(OS_LDFLAGS)

