## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard
AM_CPPFLAGS=$(GLUT_CFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/dicom -I$(top_srcdir)/unix
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mris_show
mris_show_SOURCES=mris_show.c
mris_show_LDADD= $(addprefix $(top_srcdir)/,$(LIBS_MGH)) $(GLUT_LIBS)
mris_show_LDFLAGS=
