## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard
AM_CPPFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom $(GLUT_CFLAGS) $(tcl_includes)
AM_LDFLAGS=

bin_PROGRAMS = tkmedit
tkmedit_SOURCES=tkmedit.c tkmDisplayArea.c tkmFunctionalVolume.c tkmMeditWindow.c
tkmedit_LDADD= $(LIBS_MGH) $(LIBS_TCL_OPENGL)
tkmedit_DEPENDENCIES=../utils/libutils.a
tkmedit_LDFLAGS=-L$(top_srcdir)/lib/$(OS)
