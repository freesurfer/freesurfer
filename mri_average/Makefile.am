## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=

bin_PROGRAMS = mri_average
mri_average_SOURCES=mri_average.c
mri_average_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_average_LDFLAGS=$(OS_LDFLAGS)

