## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_rigid_register
mri_rigid_register_SOURCES=mri_rigid_register.c
mri_rigid_register_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_rigid_register_LDFLAGS=$(OS_LDFLAGS)

