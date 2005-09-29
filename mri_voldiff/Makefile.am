## 
## Makefile.am 
##
AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom

bin_PROGRAMS = mri_voldiff
mri_voldiff_SOURCES=mri_voldiff.c
mri_voldiff_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_voldiff_LDFLAGS=$(OS_LDFLAGS)
