## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard
AM_CPPFLAGS=$(GLUT_CFLAGS) $(TCL_CFLAGS) $(tcl_includes) -I$(top_srcdir)/include -I$(top_srcdir)/dicom -I$(top_srcdir)/unix
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = tksurfer
tksurfer_SOURCES=tksurfer.c
tksurfer_LDADD= $(addprefix $(top_srcdir)/,$(LIBS_MGH)) $(LIBS_TCL_OPENGL)
tksurfer_LDFLAGS=
