## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_label_volume
mri_label_volume_SOURCES=mri_label_volume.c
mri_label_volume_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_label_volume_LDFLAGS=$(OS_LDFLAGS)

