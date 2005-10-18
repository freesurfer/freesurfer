## 
## Makefile.am 
##
AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom

bin_PROGRAMS = mri_morphology
mri_morphology_SOURCES=mri_morphology.c
mri_morphology_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_morphology_LDFLAGS=$(OS_LDFLAGS)
