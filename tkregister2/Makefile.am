## 
# Makefile.am - specify in addition to the standard
##

AM_CPPFLAGS=$(GLUT_CFLAGS) $(TCL_CFLAGS) $(tcl_includes) \
	-I$(top_srcdir)/include \
	-I$(top_srcdir)/dicom \
	-I$(top_srcdir)/unix \
	$(X_CFLAGS)
AM_LDFLAGS= $(LDFLAGS_TCL)

bin_PROGRAMS = tkregister2
tkregister2_SOURCES=tkregister2.c
tkregister2_LDADD= $(addprefix $(top_srcdir)/,$(LIBS_MGH)) \
	$(LIBS_TCL_OPENGL) \
	$(X_LIBS) $(GLUT_LIBS) -lGL
tkregister2_LDFLAGS=

install-exec-hook:
	mv $(bindir)/tkregister2 $(bindir)/tkregister2.bin
	cp tkregister_wrapper $(bindir)/tkregister2
	chmod a+x $(bindir)/tkregister2

EXTRA_DIST="tkregister_wrapper"

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
