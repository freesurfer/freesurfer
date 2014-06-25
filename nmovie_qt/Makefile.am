if ENABLE_QT_APPS

include $(top_srcdir)/build-aux/autotroll.mk

bin_PROGRAMS = nmovie_qt

BUILT_SOURCES = ui_MainWindow.h  \
	moc_MainWindow.cpp \
	moc_RenderWidget.cpp

nmovie_qt_SOURCES = $(BUILT_SOURCES) \
	MainWindow.cpp \
	RenderWidget.cpp \
	main.cpp

#qrc_nmovie_qt.cpp: nmovie_qt.qrc
#	$(RCC) $<  -o $@

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


nmovie_qt_CXXFLAGS = $(QT_CXXFLAGS) $(AM_CXXFLAGS)
nmovie_qt_CPPFLAGS = $(QT_CPPFLAGS) $(AM_CPPFLAGS)
nmovie_qt_LDFLAGS  = $(QT_LDFLAGS) $(OS_LDFLAGS) -framework Cocoa -framework IOKit
nmovie_qt_LDADD    = $(QT_LIBS)  \
        $(addprefix $(top_builddir)/, $(LIBS_MGH)) 

uninstall-hook:
	rm -f $(DESTDIR)$(bindir)/nmovie_qt

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


nmovie_qt_CXXFLAGS = $(QT_CXXFLAGS) $(AM_CXXFLAGS)
nmovie_qt_CPPFLAGS = $(QT_CPPFLAGS) $(AM_CPPFLAGS)
nmovie_qt_LDFLAGS  = $(QT_LDFLAGS) $(OS_LDFLAGS)
nmovie_qt_LDADD    = $(QT_LIBS)  \
        $(addprefix $(top_builddir)/, $(LIBS_MGH)) 

# put a wrapper around the bin, used to setup tcltktixblt,vtk,kww enviro vars
# and if the OS is MACOSX have a different install-hook to create a bundle
install-exec-hook:nmovie_qt
	cp $(top_builddir)/nmovie_qt/nmovie_qt $(DESTDIR)$(bindir)/nmovie_qt.bin
	echo "#!/bin/tcsh -ef" > $(DESTDIR)$(bindir)/nmovie_qt
	echo "source \$$FREESURFER_HOME/bin/tcl_setup" \
                >> $(DESTDIR)$(bindir)/nmovie_qt
	echo "source \$$FREESURFER_HOME/bin/vtk_setup" \
                >> $(DESTDIR)$(bindir)/nmovie_qt
	echo "source \$$FREESURFER_HOME/bin/qt_setup" \
                >> $(DESTDIR)$(bindir)/nmovie_qt
	echo "nmovie_qt.bin \$$argv" >> $(DESTDIR)$(bindir)/nmovie_qt
	chmod a+x $(DESTDIR)$(bindir)/nmovie_qt
	if [[ -f /usr/pubsw/bin/chrpath ]]; then \
    /usr/pubsw/bin/chrpath -d $(DESTDIR)$(bindir)/nmovie_qt.bin; fi

uninstall-hook:
	rm -f $(DESTDIR)$(bindir)/nmovie_qt.bin
endif
endif

clean-local:
	rm -f $(BUILT_SOURCES) 

EXTRA_DIST=MainWindow.ui MainWindow.h nmovie.pro
