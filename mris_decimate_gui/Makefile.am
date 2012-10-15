##
## Makefile.am - build with the GUI
##

AM_CXXFLAGS=\
	-fno-strict-aliasing \
	-Wno-deprecated \
	-I$(top_srcdir)/include \
	$(TCL_CFLAGS) \
	$(GL_CFLAGS) \
	$(VTK_CXXFLAGS) \
	-I$(top_srcdir)/vtkfsio \
	-I$(top_srcdir)/vtkutils \
	-I$(top_srcdir)/mris_decimate_gui/res \
	-I$(top_srcdir)/mris_decimate \
	$(WXWIDGETS_CXXFLAGS) $(GTS_CFLAGS)

if ENABLE_WXWIDGETS_APPS
if HAVE_GTS_LIBS
if HAVE_VTK_LIBS
    bin_PROGRAMS = mris_decimate_gui
    mris_decimate_gui_SOURCES=\
        ../mris_decimate/mris_decimate.cpp \
        ../mris_decimate/mris_decimate.h \
        DecimatePanel.cpp \
        DecimatePanel.h \
        MainApp.cpp \
        MainApp.h \
        MainWindow.cpp \
        MainWindow.h \
        MyCmdLineParser.cpp \
        MyCmdLineParser.h \
        RenderPanel.cpp \
        RenderPanel.h \
        wxVTKRenderWindowInteractor.cpp \
        wxVTKRenderWindowInteractor.h \
        res/MainWindowBase.cpp \
        res/MainWindowBase.h \
        res/MainWindowBase.fbp \
        res/DecimatePanelBase.cpp \
        res/DecimatePanelBase.h \
        res/DecimatePanelBase.fbp
    mris_decimate_gui_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH)) \
            $(VTK_LIBS) $(VTK_VERDICT_LIB) -lvtkGraphics \
            $(VTK_EXPAT_LIB) $(VTK_FREETYPE_LIB) \
            $(VTK_TIFF_LIB) $(VTK_JPEG_LIB) \
            $(VTK_METAIO_LIB) $(VTK_PNG_LIB) $(VTK_Z_LIB) \
            $(VTK_SQLITE_LIB) \
            -lvtkexoIIc -lvtkNetCDF $(VTK_NETCDFCXX_LIB) \
            -lvtkVolumeRendering -lvtkRendering -lvtkftgl \
            -lvtkHybrid $(VTK_HDF5_LIB) \
						-lvtkWidgets \
            -lvtkGraphics \
            -lvtkCommon  \
            -lvtksys \
            -lvtkImaging \
            -lvtkFiltering \
            -lvtkRendering \
            -lvtkIO $(VTK_LSDYNA_LIB) -lvtkDICOMParser \
            $(WXWIDGETS_LIBS) $(WXWIDGETS_GL_LIBS) $(GTS_LIBS)
    mris_decimate_gui_LDFLAGS=$(OS_LDFLAGS)

# put a wrapper around the bin, used to setup vtk enviro vars

# and if the OS is MACOSX have a different install-hook to create a bundle
if HAVE_MAC_OSX
install-exec-hook:mris_decimate_gui
	mkdir -p $(DESTDIR)$(bindir)/MrisDecimate.app
	mkdir -p $(DESTDIR)$(bindir)/MrisDecimate.app/Contents
	mkdir -p $(DESTDIR)$(bindir)/MrisDecimate.app/Contents/MacOS
	mkdir -p $(DESTDIR)$(bindir)/MrisDecimate.app/Contents/Resources
	mkdir -p $(DESTDIR)$(bindir)/MrisDecimate.app/Contents/Resources/English.lproj
	cp $(top_builddir)/mris_decimate_gui/Info.plist $(DESTDIR)$(bindir)/MrisDecimate.app/Contents
	echo -n 'APPL????' > $(DESTDIR)$(bindir)/MrisDecimate.app/Contents/PkgInfo
	cp $(top_builddir)/mris_decimate_gui/mris_decimate_gui $(DESTDIR)$(bindir)/MrisDecimate.app/Contents/MacOS/mris_decimate_gui.bin
	echo "#!/bin/tcsh -ef" > $(DESTDIR)$(bindir)/mris_decimate_gui
	echo "source \$$FREESURFER_HOME/bin/tcl_setup" \
		>> $(DESTDIR)$(bindir)/mris_decimate_gui
	echo "source \$$FREESURFER_HOME/bin/vtk_setup" \
		>> $(DESTDIR)$(bindir)/mris_decmris_deci
	echo "exec $(DESTDIR)$(bindir)/MrisDecimate.app/Contents/MacOS/mris_decimate_gui.bin \$$*" >> $(DESTDIR)$(bindir)/mris_decimate_gui
	chmod a+x $(DESTDIR)$(bindir)/mris_decimate_gui
	if [[ -f /usr/pubsw/bin/chrpath ]]; then \
	/usr/pubsw/bin/chrpath -d $(DESTDIR)$(bindir)/MrisDecimate.app/Contents/MacOS/mris_decimate_gui.bin; fi

uninstall-hook:
	rm -f $(DESTDIR)$(bindir)/mris_decimate_gui
	rm -rf $(DESTDIR)$(bindir)/MrisDecimatey.app
else
install-exec-hook:mris_decimate_gui
	cp $(top_builddir)/mris_decimate_gui/mris_decimate_gui $(DESTDIR)$(bindir)/mris_decimate_gui.bin
	echo "#!/bin/tcsh -ef" > $(DESTDIR)$(bindir)/mris_decimate_gui
	echo "source \$$FREESURFER_HOME/bin/vtk_setup" \
		>> $(DESTDIR)$(bindir)/mris_decimate_gui
	echo "if (-e /bin/pwd) then" >> $(DESTDIR)$(bindir)/mris_decimate_gui
	echo "  setenv HOME \`cd \$$HOME; /bin/pwd\`" \
		>> $(DESTDIR)$(bindir)/mris_decimate_gui
	echo "  setenv SUBJECTS_DIR \`cd \$$SUBJECTS_DIR; /bin/pwd\`" \
		>> $(DESTDIR)$(bindir)/mris_decimate_gui
	echo "endif" >> $(DESTDIR)$(bindir)/mris_decimate_gui
	echo "mris_decimate_gui.bin \$$argv" >> $(DESTDIR)$(bindir)/mris_decimate_gui
	chmod a+x $(DESTDIR)$(bindir)/mris_decimate_gui
	if [[ -f /usr/pubsw/bin/chrpath ]]; then \
	/usr/pubsw/bin/chrpath -d $(DESTDIR)$(bindir)/mris_decimate_gui.bin; fi

uninstall-hook:
	rm -f $(DESTDIR)$(bindir)/mris_decimate_gui.bin
endif

endif
endif
endif

EXTRA_DIST=Info.plist

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
