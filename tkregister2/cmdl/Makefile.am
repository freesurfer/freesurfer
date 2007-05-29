## 
# Makefile.am
##

# build without Tcl/Tk and OpenGL (command-line processing only)
AM_CPPFLAGS=-DNO_GUI -I$(top_srcdir)/include
AM_LDFLAGS=
bin_PROGRAMS = tkregister2_cmdl
tkregister2_cmdl_SOURCES=../tkregister2.c
tkregister2_cmdl_LDADD= $(addprefix $(top_srcdir)/,$(LIBS_MGH))

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
