## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_relabel_hypointensities
mri_relabel_hypointensities_SOURCES=mri_relabel_hypointensities.c
mri_relabel_hypointensities_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_relabel_hypointensities_LDFLAGS=$(OS_LDFLAGS)

