## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=

bin_PROGRAMS = mri_build_priors
mri_build_priors_SOURCES=mri_build_priors.c
mri_build_priors_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_build_priors_LDFLAGS=$(OS_LDFLAGS)

