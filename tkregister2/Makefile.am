## 
# Makefile.am
##

# also build command-line only version (no gui)
SUBDIRS=cmdl

if ENABLE_TCLTK_APPS

AM_CPPFLAGS=$(TCL_CFLAGS) $(tcl_includes) \
	-I$(top_srcdir)/include \
	-I$(top_srcdir)/unix \
	-I$(top_srcdir)/glut \
	$(GL_CFLAGS)
AM_LDFLAGS= $(LDFLAGS_TCL)
bin_PROGRAMS = tkregister2
tkregister2_SOURCES=tkregister2.c
tkregister2_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH)) \
	$(LIBS_TCL_OPENGL) \
	$(X_LIBS) -lGL \
	$(top_builddir)/glut/libglut.a

else

bin_PROGRAMS = 

endif

if ENABLE_TCLTK_APPS
# put a wrapper around the bin, used to setup tcltktixblt enviro vars
# and install in 'tktools' dir instead of 'bin' (to allow easy deletion
# of all cortech-licensed software)
install-exec-hook:
	mv -v $(DESTDIR)$(bindir)/tkregister2_cmdl $(DESTDIR)$(prefix)/tktools/
	rm -vf $(DESTDIR)$(bindir)/tkregister*
	cp $(top_builddir)/tkregister2/tkregister2 \
		$(DESTDIR)$(prefix)/tktools/tkregister2.bin
	echo "#!/bin/tcsh -ef" > $(DESTDIR)$(prefix)/tktools/tkregister2
	echo "source \$$FREESURFER_HOME/bin/tcl_setup" \
		>> $(DESTDIR)$(prefix)/tktools/tkregister2
	echo "exec \$$FREESURFER_HOME/tktools/tkregister2.bin \$$argv" \
		>> $(DESTDIR)$(prefix)/tktools/tkregister2
	cp $(DESTDIR)$(prefix)/tktools/tkregister2 \
		$(DESTDIR)$(prefix)/tktools/tkregister
	chmod a+x $(DESTDIR)$(prefix)/tktools/tkregister
	chmod a+x $(DESTDIR)$(prefix)/tktools/tkregister2
	cp $(DESTDIR)$(prefix)/tktools/tkregister $(DESTDIR)$(bindir)/
	cp $(DESTDIR)$(prefix)/tktools/tkregister2 $(DESTDIR)$(bindir)/
	if [[ -f /usr/pubsw/bin/chrpath ]]; then \
	/usr/pubsw/bin/chrpath -d $(DESTDIR)$(prefix)/tktools/tkregister2.bin; fi

uninstall-hook:
	rm -f $(DESTDIR)$(prefix)/tktools/tkregister*
	rm -f $(DESTDIR)$(bindir)/tkregister*

else
all:
	@echo ""
	@echo "Excluding Tcl/Tk and OpenGL elements from this tkregister2 build."
	@echo "Run configure with either --enable-tcltk-apps or --with-tcl-dir=DIR"
	@echo ""
endif

# Copy over precompiled 32-bit versions until we get 64-bit builds working.
if INSTALL_MAC32_GUIS
install-exec-hook:
	mv -v $(DESTDIR)$(bindir)/tkregister2_cmdl $(DESTDIR)$(prefix)/tktools/
	rm -vf $(DESTDIR)$(bindir)/tkregister*
	cp $(top_builddir)/distribution/Mac_32b/tkregister2_Mac32.bin \
		$(DESTDIR)$(prefix)/tktools/tkregister2.bin
	echo "#!/bin/tcsh -ef" > $(DESTDIR)$(prefix)/tktools/tkregister2
	echo "source \$$FREESURFER_HOME/bin/tcl_setup" \
		>> $(DESTDIR)$(prefix)/tktools/tkregister2
	echo "exec \$$FREESURFER_HOME/tktools/tkregister2.bin \$$argv" \
		>> $(DESTDIR)$(prefix)/tktools/tkregister2
	cp $(DESTDIR)$(prefix)/tktools/tkregister2 \
		$(DESTDIR)$(prefix)/tktools/tkregister
	chmod a+x $(DESTDIR)$(prefix)/tktools/tkregister
	chmod a+x $(DESTDIR)$(prefix)/tktools/tkregister2
	cp $(DESTDIR)$(prefix)/tktools/tkregister $(DESTDIR)$(bindir)/
	cp $(DESTDIR)$(prefix)/tktools/tkregister2 $(DESTDIR)$(bindir)/
	if [[ -f /usr/pubsw/bin/chrpath ]]; then \
	/usr/pubsw/bin/chrpath -d $(DESTDIR)$(prefix)/tktools/tkregister2.bin; fi

uninstall-hook:
	rm -f $(DESTDIR)$(prefix)/tktools/tkregister*
	rm -f $(DESTDIR)$(bindir)/tkregister*
endif



# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
