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
	$(TCL_CFLAGS) $(TIFFCFLAGS) \
	$(tcl_includes) $(GL_CFLAGS) \
	-I$(top_srcdir)/include \
	-I$(top_srcdir)/glut
AM_LDFLAGS=$(LDFLAGS_TCL) $(TIFF_LIBS)

if ENABLE_TCLTK_APPS
bin_PROGRAMS=tkmedit
tkmedit_SOURCES=tkmedit.c \
	tkmDisplayArea.c \
	tkmFunctionalVolume.c \
	tkmMeditWindow.c \
	tkmedit.h \
	tkmDisplayArea.h \
	tkmFunctionalVolume.h \
	tkmMeditWindow.h \
	$(top_srcdir)/fsgdf/fsgdf_wrap.c
tkmedit_LDADD= $(top_builddir)/fsgdf/libfsgdf.a \
	$(addprefix $(top_builddir)/,$(LIBS_MGH)) \
	$(LIBS_TCL_OPENGL) \
	$(X_LIBS) $(LIBTIFF) \
	$(top_builddir)/glut/libglut.a
tkmedit_LDFLAGS=$(LDFLAGS_TCL)
else
bin_PROGRAMS=
endif

if ENABLE_TCLTK_APPS
# put a wrapper around the bin, used to setup tcltktixblt enviro vars
# and install in 'tktools' dir instead of 'bin' (to allow easy deletion
# of all cortech-licensed software)
install-exec-hook:
	rm -vf $(DESTDIR)$(bindir)/tkmedit
	cp $(top_builddir)/tkmedit/tkmedit $(DESTDIR)$(prefix)/tktools/tkmedit.bin
	echo "#!/bin/tcsh -ef" > $(DESTDIR)$(prefix)/tktools/tkmedit
	echo "source \$$FREESURFER_HOME/bin/tcl_setup" \
		>> $(DESTDIR)$(prefix)/tktools/tkmedit
	echo "exec \$$FREESURFER_HOME/tktools/tkmedit.bin \$$argv" \
		>> $(DESTDIR)$(prefix)/tktools/tkmedit
	chmod a+x $(DESTDIR)$(prefix)/tktools/tkmedit
	cp $(DESTDIR)$(prefix)/tktools/tkmedit $(DESTDIR)$(bindir)/
	if [[ -f /usr/pubsw/bin/chrpath ]]; then \
	/usr/pubsw/bin/chrpath -d $(DESTDIR)$(prefix)/tktools/tkmedit.bin; fi

uninstall-hook:
	rm -f $(DESTDIR)$(prefix)/tktools/tkmedit*

else
all:
	@echo "Building Tcl/Tk apps not enabled."
	@echo "Run configure with either --enable-tcltk-apps or --with-tcl-dir=DIR"
install-exec-hook:
endif

# Copy over precompiled 32-bit versions until we get 64-bit builds working.
if INSTALL_MAC32_GUIS
install-exec-hook:
	rm -vf $(DESTDIR)$(bindir)/tkmedit
	cp $(top_builddir)/distribution/Mac_32b/tkmedit_Mac32.bin $(DESTDIR)$(prefix)/tktools/tkmedit.bin
	echo "#!/bin/tcsh -ef" > $(DESTDIR)$(prefix)/tktools/tkmedit
	echo "source \$$FREESURFER_HOME/bin/tcl_setup" \
		>> $(DESTDIR)$(prefix)/tktools/tkmedit
	echo "exec \$$FREESURFER_HOME/tktools/tkmedit.bin \$$argv" \
		>> $(DESTDIR)$(prefix)/tktools/tkmedit
	chmod a+x $(DESTDIR)$(prefix)/tktools/tkmedit
	cp $(DESTDIR)$(prefix)/tktools/tkmedit $(DESTDIR)$(bindir)/
	if [[ -f /usr/pubsw/bin/chrpath ]]; then \
	/usr/pubsw/bin/chrpath -d $(DESTDIR)$(prefix)/tktools/tkmedit.bin; fi

uninstall-hook:
	rm -f $(DESTDIR)$(prefix)/tktools/tkmedit*
endif

EXTRA_DIST=development_notes.txt

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
