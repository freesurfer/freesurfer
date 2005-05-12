## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom

bin_PROGRAMS = mri_vol2vol
mri_vol2vol_SOURCES=mri_vol2vol.c
mri_vol2vol_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_vol2vol_LDFLAGS=$(OS_LDFLAGS)

