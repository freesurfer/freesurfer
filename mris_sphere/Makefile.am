## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mris_sphere
mris_sphere_SOURCES=mris_sphere.c
mris_sphere_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_sphere_LDFLAGS=$(OS_LDFLAGS)

