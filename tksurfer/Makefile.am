## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard
AM_CPPFLAGS=$(GLUT_CFLAGS) $(TCL_CFLAGS) $(tcl_includes) -I$(top_srcdir)/include -I$(top_srcdir)/dicom -I$(top_srcdir)/unix
AM_LDFLAGS=$(LDFLAGS_TCL)

bin_PROGRAMS = tksurfer
tksurfer_SOURCES=tksurfer.c $(top_srcdir)/fsgdf/fsgdf_wrap.c
tksurfer_LDADD= $(addprefix $(top_srcdir)/,$(LIBS_MGH)) $(top_srcdir)/fsgdf/libfsgdf.a $(LIBS_TCL_OPENGL) -lpthread
tksurfer_LDFLAGS=
