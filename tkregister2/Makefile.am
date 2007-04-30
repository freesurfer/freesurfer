## 
# Makefile.am
##

if ENABLE_TCLTK_APPS

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

else

# build without Tcl/Tk and OpenGL (command-line processing only)
AM_CPPFLAGS=-DNO_GUI
	-I$(top_srcdir)/include \
	-I$(top_srcdir)/dicom \
	-I$(top_srcdir)/unix
AM_LDFLAGS=
bin_PROGRAMS = tkregister2
tkregister2_SOURCES=tkregister2.c
tkregister2_LDADD= $(addprefix $(top_srcdir)/,$(LIBS_MGH))
tkregister2_LDFLAGS=

endif

if ENABLE_TCLTK_APPS
install-exec-hook:
	cp $(srcdir)/tkregister2 $(DESTDIR)$(bindir)/tkregister2.bin
	echo "#!/bin/tcsh -ef" > $(DESTDIR)$(bindir)/tkregister2
	echo "source \$$FREESURFER_HOME/bin/tcl_setup" \
		>> $(DESTDIR)$(bindir)/tkregister2
	echo "tkregister2.bin \$$argv" >> $(DESTDIR)$(bindir)/tkregister2
	chmod a+x $(DESTDIR)$(bindir)/tkregister2
else
all:
	@echo "Excluding Tcl/Tk and OpenGL elements from this tkregister2 build."
  @echo "Run configure with either --enable-tcltk-apps or --with-tcl-dir=DIR"
endif

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
