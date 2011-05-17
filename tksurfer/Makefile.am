## 
# Makefile.am
##

if USE_LOCAL_TIFF
TIFFCFLAGS=-I$(top_srcdir)/tiff
LIBTIFF=$(top_builddir)/tiff/libtiff.a
else
TIFFCFLAGS=$(TIFF_CFLAGS)
LIBTIFF=$(LIB_TIFF)
endif

AM_CPPFLAGS=\
	$(TCL_CFLAGS) $(GL_CFLAGS) $(TIFFCFLAGS) \
	-I$(top_srcdir)/include \
	-I$(top_srcdir)/dicom \
	-I$(top_srcdir)/glut \
	-I$(top_srcdir)/unix

AM_LDFLAGS=$(LDFLAGS_TCL) $(TIFF_LIBS)

if ENABLE_TCLTK_APPS
bin_PROGRAMS=tksurfer
tksurfer_SOURCES=tksurfer.c
tksurfer_LDADD= $(top_builddir)/fsgdf/libfsgdf.a \
  $(addprefix $(top_builddir)/, $(LIBS_MGH)) \
	$(LIBS_TCL_OPENGL) \
	$(X_LIBS) -lGL $(LIBTIFF) \
	$(top_builddir)/glut/libglut.a
tksurfer_LDFLAGS=$(LDFLAGS_TCL)
else
bin_PROGRAMS=
endif

if ENABLE_TCLTK_APPS
# put a wrapper around the bin, used to setup tcltktixblt enviro vars
# and install in 'tktools' dir instead of 'bin' (to allow easy deletion
# of all cortech-licensed software)
install-exec-hook:
	rm -fv $(DESTDIR)$(bindir)/tksurfer
	cp $(top_builddir)/tksurfer/tksurfer $(DESTDIR)$(prefix)/tktools/tksurfer.bin
	echo "#!/bin/tcsh -ef" > $(DESTDIR)$(prefix)/tktools/tksurfer
	echo "source \$$FREESURFER_HOME/bin/tcl_setup" \
		>> $(DESTDIR)$(prefix)/tktools/tksurfer
	echo "exec \$$FREESURFER_HOME/tktools/tksurfer.bin \$$argv" \
		>> $(DESTDIR)$(prefix)/tktools/tksurfer
	chmod a+x $(DESTDIR)$(prefix)/tktools/tksurfer
	cp $(DESTDIR)$(prefix)/tktools/tksurfer $(DESTDIR)$(bindir)/
	if [[ -f /usr/pubsw/bin/chrpath ]]; then \
	/usr/pubsw/bin/chrpath -d $(DESTDIR)$(prefix)/tktools/tksurfer.bin; fi

uninstall-hook:
	rm -f $(DESTDIR)$(prefix)/tktools/tksurfer*

else
all:
	@echo "Building Tcl/Tk apps not enabled."
	@echo "Run configure with either --enable-tcltk-apps or --with-tcl-dir=DIR"
install-exec-hook:
endif

EXTRA_DIST=development_notes.txt

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
