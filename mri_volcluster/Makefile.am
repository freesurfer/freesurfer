## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_volcluster
mri_volcluster_SOURCES=mri_volcluster.c
mri_volcluster_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_volcluster_LDFLAGS=$(OS_LDFLAGS)

