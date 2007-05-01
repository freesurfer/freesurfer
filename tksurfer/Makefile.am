## 
# Makefile.am
##

AM_CPPFLAGS=$(GLUT_CFLAGS) \
	$(TCL_CFLAGS) \
	$(X_CFLAGS) \
	-I$(top_srcdir)/include \
	-I$(top_srcdir)/dicom \
	-I$(top_srcdir)/unix

AM_LDFLAGS=$(LDFLAGS_TCL)

if ENABLE_TCLTK_APPS
bin_PROGRAMS=tksurfer
tksurfer_SOURCES=tksurfer.c $(top_srcdir)/fsgdf/fsgdf_wrap.c
tksurfer_LDADD= 	$(top_srcdir)/fsgdf/libfsgdf.a \
  $(addprefix $(top_srcdir)/,$(LIBS_MGH)) \
	$(LIBS_TCL_OPENGL) \
	$(X_LIBS) $(GLUT_LIBS) -lGL
tksurfer_LDFLAGS=$(LDFLAGS_TCL)
else
bin_PROGRAMS=
endif

if ENABLE_TCLTK_APPS
install-exec-hook:
	cp $(srcdir)/tksurfer $(DESTDIR)$(bindir)/tksurfer.bin
	echo "#!/bin/tcsh -ef" > $(DESTDIR)$(bindir)/tksurfer
	echo "source \$$FREESURFER_HOME/bin/tcl_setup" \
		>> $(DESTDIR)$(bindir)/tksurfer
	echo "tksurfer.bin \$$argv" >> $(DESTDIR)$(bindir)/tksurfer
	chmod a+x $(DESTDIR)$(bindir)/tksurfer
else
all:
	@echo "Building Tcl/Tk apps not enabled."
	@echo "Run configure with either --enable-tcltk-apps or --with-tcl-dir=DIR"
install-exec-hook:
endif

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
