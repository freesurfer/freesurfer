#
# Makefile.am for qdec
#

SUBDIRS=images

AM_CPPFLAGS=-I$(top_srcdir)/include

AM_LDFLAGS=

FORMS = mainform.ui

UICDECLS = .ui/mainform.h mainform.ui.h

UICIMPLS = .ui/mainform.cpp

SRCMOC   = .moc/moc_mainform.cpp

MOC      = $(QTDIR)/bin/moc
UIC      = $(QTDIR)/bin/uic

#don't build Qt-based apps
bin_PROGRAMS=
#bin_PROGRAMS=qdec ${QDEC_BUNDLE_TARGET}
EXTRA_PROGRAMS=

qdec_SOURCES=\
	Utilities.cpp \
	Utilities.h \
	$(FORMS) \
	$(UICDECLS) \
	$(UICIMPLS) \
	$(SRCMOC) \
	$(OBJMOC) \
	.ui/image_collection.cpp \
	main.cpp 
qdec_LDADD=$(addprefix $(top_srcdir)/,$(LIBS_MGH)) \
	-L$(QTDIR)/lib -lqt-mt
qdec_LDFLAGS=
qdec_CXXFLAGS=-DQT_SHARED \
	-DQT_NO_DEBUG \
	-DQT_THREAD_SUPPORT \
	-I$(QTDIR)/include \
	-I.ui/ 

# All the Qt stuff
.ui/mainform.h: mainform.ui
	$(UIC) mainform.ui -o $@
.ui/mainform.cpp: \
	.ui/mainform.h \
	mainform.ui
	$(UIC) mainform.ui -i mainform.h -o $@
mainform.$(OBJEXT): .ui/mainform.cpp 
	mainform.ui.h .ui/mainform.h

moc_mainform.$(OBJEXT): .moc/moc_mainform.cpp
.moc/moc_mainform.cpp: .ui/mainform.h
	$(MOC) .ui/mainform.h -o $@

image_collection.$(OBJEXT): .ui/image_collection.cpp

.ui/image_collection.cpp: \
		images/fileopen \
		images/filesave \
		images/print \
		images/undo \
		images/redo \
		images/editcut \
		images/editcopy \
		images/editpaste \
		images/searchfind
	$(UIC)  -embed qdec \
	images/filenew \
	images/fileopen \
	images/filesave \
	images/print \
	images/undo \
	images/redo \
	images/editcut \
	images/editcopy \
	images/editpaste \
	images/searchfind \
	-o $@

qdec_bundle: qdec
	@echo Making OSX bundle...
	test -d qdec.app/Contents/MacOS/ || mkdir -p qdec.app/Contents/MacOS/
	rm -f qdec.app/Contents/MacOS/../PkgInfo
	@echo "APPL????" >qdec.app/Contents/MacOS/../PkgInfo
	test -d qdec.app/Contents/MacOS/ || mkdir -p qdec.app/Contents/MacOS/
	rm -f qdec.app/Contents/MacOS/../Info.plist
	sed -e "s,@ICON@,application.icns,g" -e "s,@EXECUTABLE@,qdec,g" "$(QTDIR)/mkspecs/default/Info.plist.app" >"qdec.app/Contents/MacOS/../Info.plist"
	test -d qdec.app/Contents/MacOS/ || mkdir -p qdec.app/Contents/MacOS/
	rm -f qdec.app/Contents/MacOS/qdec
	cp qdec qdec.app/Contents/MacOS/
	test -d qdec.app/Contents/Frameworks || mkdir -p qdec.app/Contents/Frameworks
	cp $(QTDIR)/lib/libqt-mt.3.dylib qdec.app/Contents/Frameworks
	install_name_tool -id @executable_path/../Frameworks/libqt-mt.3.dylib qdec.app/Contents/Frameworks/libqt-mt.3.dylib
	install_name_tool -change libqt-mt.3.dylib @executable_path/../Frameworks/libqt-mt.3.dylib qdec.app/Contents/MacOS/qdec
	cp -r qdec.app $(DESTDIR)$(bindir)
	echo "$(FREESURFER_HOME)/bin/qdec.app/Contents/MacOS/qdec" > $(DESTDIR)$(bindir)/qdec
	chmod a+x $(DESTDIR)$(bindir)/qdec

.PHONY: qdec_bundle

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
