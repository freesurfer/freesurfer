## 
## Makefile.am 
##
#######################################################################
AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom

bin_PROGRAMS = mri_cht2p
mri_cht2p_SOURCES=mri_cht2p.c
mri_cht2p_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH)) -lgsl -lgslcblas
mri_cht2p_LDFLAGS=$(OS_LDFLAGS)




