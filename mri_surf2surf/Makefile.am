## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mri_surf2surf
mri_surf2surf_SOURCES=mri_surf2surf.c
mri_surf2surf_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_surf2surf_LDFLAGS=$(OS_LDFLAGS)

