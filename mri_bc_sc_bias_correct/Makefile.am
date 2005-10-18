## 
## Makefile.am 
##
AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom

bin_PROGRAMS = mri_bc_sc_bias_correct
mri_bc_sc_bias_correct_SOURCES=mri_bc_sc_bias_correct.c
mri_bc_sc_bias_correct_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_bc_sc_bias_correct_LDFLAGS=$(OS_LDFLAGS)
