## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_vol2surf
mri_vol2surf_SOURCES=mri_vol2surf.c
mri_vol2surf_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_vol2surf_LDFLAGS=$(OS_LDFLAGS)

