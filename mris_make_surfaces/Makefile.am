## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mris_make_surfaces mris_refine_surfaces
mris_make_surfaces_SOURCES=mris_make_surfaces.c
mris_make_surfaces_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_make_surfaces_LDFLAGS= $(OS_LDFLAGS)

mris_refine_surfaces_SOURCES=mris_refine_surfaces.c
mris_refine_surfaces_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_refine_surfaces_LDFLAGS= $(OS_LDFLAGS)

