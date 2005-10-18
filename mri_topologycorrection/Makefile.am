## 
## Makefile.am 
##
AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom

bin_PROGRAMS = mri_topologycorrection
mri_topologycorrection_SOURCES=mri_topologycorrection.c
mri_topologycorrection_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_topologycorrection_LDFLAGS=$(OS_LDFLAGS)
