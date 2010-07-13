##
## Makefile.am 
##

SUBDIRS=gui

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include $(VTK_CXXFLAGS) -Wno-deprecated

if HAVE_VTK_LIBS
    bin_PROGRAMS = mris_decimate
    mris_decimate_SOURCES=mris_decimate.cpp mris_decimate.h
    mris_decimate_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH)) \
				$(VTK_LIBS) $(VTK_VERDICT_LIB) -lvtkHybrid -lvtkWidgets \
        -lvtkGraphics \
        -lvtkCommon  \
        -lvtkImaging \
        -lvtkFiltering
    mris_decimate_LDFLAGS=$(OS_LDFLAGS)

# put a wrapper around the bin, used to setup vtk enviro vars
install-exec-hook:mris_decimate mris_decimate_gui
	cp $(top_builddir)/mris_decimate/mris_decimate $(DESTDIR)$(bindir)/mris_decimate.bin
	echo "#!/bin/tcsh -ef" > $(DESTDIR)$(bindir)/mris_decimate
	echo "source \$$FREESURFER_HOME/bin/vtk_setup" \
		>> $(DESTDIR)$(bindir)/mris_decimate
	echo "if (-e /bin/pwd) then" >> $(DESTDIR)$(bindir)/mris_decimate
	echo "  setenv HOME \`cd \$$HOME; /bin/pwd\`" \
		>> $(DESTDIR)$(bindir)/mris_decimate
	echo "  setenv SUBJECTS_DIR \`cd \$$SUBJECTS_DIR; /bin/pwd\`" \
		>> $(DESTDIR)$(bindir)/mris_decimate
	echo "endif" >> $(DESTDIR)$(bindir)/mris_decimate
	echo "mris_decimate.bin \$$argv" >> $(DESTDIR)$(bindir)/mris_decimate
	chmod a+x $(DESTDIR)$(bindir)/mris_decimate
	if [[ -f /usr/pubsw/bin/chrpath ]]; then \
	/usr/pubsw/bin/chrpath -d $(DESTDIR)$(bindir)/mris_decimate.bin; fi

uninstall-hook:
	rm -f $(DESTDIR)$(bindir)/mris_decimate.bin

endif

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
