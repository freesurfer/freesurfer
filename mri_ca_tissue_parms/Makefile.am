## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=

bin_PROGRAMS = mri_ca_tissue_parms
mri_ca_tissue_parms_SOURCES=mri_ca_tissue_parms.c
mri_ca_tissue_parms_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_ca_tissue_parms_LDFLAGS=$(OS_LDFLAGS)

