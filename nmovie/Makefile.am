## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard
AM_CPPFLAGS=$(GLUT_CFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/dicom -I$(top_srcdir)/unix
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = nmovie
nmovie_SOURCES=nmovie.c
nmovie_LDADD= $(addprefix $(top_srcdir)/,$(LIBS_MGH)) -lXt -lXext -lXaw $(GLUT_LIBS)  
nmovie_LDFLAGS=
