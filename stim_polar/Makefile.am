## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=$(GLUT_CFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = stim_polar
stim_polar_SOURCES=stim_polar.c $(top_srcdir)/utils/machine.c
stim_polar_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH)) $(GLUT_LIBS)
stim_polar_LDFLAGS=

