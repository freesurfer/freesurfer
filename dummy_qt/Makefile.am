if ENABLE_QT_APPS

include $(top_srcdir)/build-aux/autotroll.mk

bin_PROGRAMS = dummy_qt

BUILT_SOURCES = ui_MainWindow.h  \
	moc_MainWindow.cpp \
	qrc_dummy_qt.cpp

dummy_qt_SOURCES = $(BUILT_SOURCES) \
	MainWindow.cpp \
	main.cpp \
	dummy_qt.qrc

qrc_dummy_qt.cpp: dummy_qt.qrc
	$(RCC) $<  -o $@

if HAVE_MAC_OSX
AM_CXXFLAGS=\
	-fno-strict-aliasing \
	-Wno-deprecated \
	-Wno-unused \
	-Wno-uninitialized \
	-Wno-return-type \
	-Wno-reorder \
	-Wno-sign-compare \
	-I$(top_srcdir)/include \
	$(GL_CFLAGS) 


dummy_qt_CXXFLAGS = $(QT_CXXFLAGS) $(AM_CXXFLAGS)
dummy_qt_CPPFLAGS = $(QT_CPPFLAGS) $(AM_CPPFLAGS)
dummy_qt_LDFLAGS  = $(QT_LDFLAGS) $(OS_LDFLAGS) -framework Cocoa -framework IOKit
dummy_qt_LDADD    = $(QT_LIBS)  \
        $(addprefix $(top_builddir)/, $(LIBS_MGH)) 

# need to create a bundle for Macs
install-exec-hook:dummy_qt
	mkdir -p $(DESTDIR)$(bindir)/DummyQt.app
	mkdir -p $(DESTDIR)$(bindir)/DummyQt.app/Contents
	mkdir -p $(DESTDIR)$(bindir)/DummyQt.app/Contents/MacOS
	mkdir -p $(DESTDIR)$(bindir)/DummyQt.app/Contents/Resources
	mkdir -p $(DESTDIR)$(bindir)/DummyQt.app/Contents/Resources/English.lproj
#       mkdir -p $(DESTDIR)$(bindir)/DummyQt.app/Contents/Frameworks
	cp $(top_builddir)/dummy_qt/Info.plist $(DESTDIR)$(bindir)/DummyQt.app/Contents
	echo -n 'APPL????' > $(DESTDIR)$(bindir)/DummyQt.app/Contents/PkgInfo
	cp $(top_builddir)/dummy_qt/dummy_qt $(DESTDIR)$(bindir)/DummyQt.app/Contents/MacOS/dummy_qt.bin
	echo "#!/bin/tcsh -ef" > $(DESTDIR)$(bindir)/dummy_qt
#       echo "source \$$FREESURFER_HOME/bin/tcl_setup" \
                >> $(DESTDIR)$(bindir)/dummy_qt
	echo "source \$$FREESURFER_HOME/bin/qt_setup" \
                >> $(DESTDIR)$(bindir)/dummy_qt
	echo "exec $(DESTDIR)$(bindir)/DummyQt.app/Contents/MacOS/dummy_qt.bin \$$*" >> $(DESTDIR)$(bindir)/dummy_qt
	chmod a+x $(DESTDIR)$(bindir)/dummy_qt
#       ./mac_package_dylibs.sh $(WXWIDGETS_DIR)lib $(DESTDIR)$(bindir)/Freeview.app/Contents/Frameworks $(DESTDIR)$(bindir)/Freeview.app/Contents/MacOS/freeview.bin
	if [[ -f /usr/pubsw/bin/chrpath ]]; then \
	/usr/pubsw/bin/chrpath -d $(DESTDIR)$(bindir)/DummyQt.App/Contents/MacOS/dummy_qt.bin; fi

uninstall-hook:
	rm -f $(DESTDIR)$(bindir)/dummy_qt
	rm -rf $(DESTDIR)$(bindir)/DummyQt.app

else

AM_CXXFLAGS=\
	-fno-strict-aliasing \
	-Wno-deprecated \
	-Wno-unused \
	-Wno-uninitialized \
	-Wno-return-type \
	-Wno-reorder \
	-Wno-sign-compare \
	-I$(top_srcdir)/include \
	$(GL_CFLAGS) 


dummy_qt_CXXFLAGS = $(QT_CXXFLAGS) $(AM_CXXFLAGS)
dummy_qt_CPPFLAGS = $(QT_CPPFLAGS) $(AM_CPPFLAGS)
dummy_qt_LDFLAGS  = $(QT_LDFLAGS) $(OS_LDFLAGS) 
dummy_qt_LDADD    = $(QT_LIBS)  \
        $(addprefix $(top_builddir)/, $(LIBS_MGH)) 

# put a wrapper around the bin, used to setup tcltktixblt,vtk,kww enviro vars
# and if the OS is MACOSX have a different install-hook to create a bundle
install-exec-hook:dummy_qt
	cp $(top_builddir)/dummy_qt/dummy_qt $(DESTDIR)$(bindir)/dummy_qt.bin
	echo "#!/bin/tcsh -ef" > $(DESTDIR)$(bindir)/dummy_qt
	echo "source \$$FREESURFER_HOME/bin/tcl_setup" \
                >> $(DESTDIR)$(bindir)/dummy_qt
	echo "source \$$FREESURFER_HOME/bin/vtk_setup" \
                >> $(DESTDIR)$(bindir)/dummy_qt
	echo "source \$$FREESURFER_HOME/bin/qt_setup" \
                >> $(DESTDIR)$(bindir)/dummy_qt
	echo "dummy_qt.bin \$$argv" >> $(DESTDIR)$(bindir)/dummy_qt
	chmod a+x $(DESTDIR)$(bindir)/dummy_qt
	if [[ -f /usr/pubsw/bin/chrpath ]]; then \
	/usr/pubsw/bin/chrpath -d $(DESTDIR)$(bindir)/dummy_qt.bin; fi

uninstall-hook:
	rm -f $(DESTDIR)$(bindir)/dummy_qt.bin
endif
endif

EXTRA_DIST=MainWindow.ui MainWindow.h dummy_qt.pro
