## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard
AM_CPPFLAGS=$(GLUT_CFLAGS) $(tcl_includes) -I$(top_srcdir)/include -I$(top_srcdir)/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = tkmedit
tkmedit_SOURCES=tkmedit.c tkmDisplayArea.c tkmFunctionalVolume.c tkmMeditWindow.c
tkmedit_LDADD= $(addprefix $(top_srcdir)/,$(LIBS_MGH)) $(LIBS_TCL_OPENGL)
tkmedit_DEPENDENCIES=$(addprefix $(top_srcdir)/,$(LIBS_MGH))

