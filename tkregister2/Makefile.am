## 
# Makefile.am - specify in addition to the standard
##

AM_CPPFLAGS=$(GLUT_CFLAGS) $(TCL_CFLAGS) $(tcl_includes) \
	-I$(top_srcdir)/include \
	-I$(top_srcdir)/dicom \
	-I$(top_srcdir)/unix \
	$(X_CFLAGS)
AM_LDFLAGS= $(LDFLAGS_TCL)

if ENABLE_TCLTK_APPS
bin_PROGRAMS = tkregister2
else
bin_PROGRAMS =
endif
tkregister2_SOURCES=tkregister2.c
tkregister2_LDADD= $(addprefix $(top_srcdir)/,$(LIBS_MGH)) \
	$(LIBS_TCL_OPENGL) \
	$(X_LIBS) $(GLUT_LIBS) -lGL
tkregister2_LDFLAGS=


# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
