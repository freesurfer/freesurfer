## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_tessellate
mri_tessellate_SOURCES=mri_tessellate.c
mri_tessellate_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_tessellate_LDFLAGS=$(OS_LDFLAGS)

