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

bin_PROGRAMS=qdec
EXTRA_PROGRAMS=

qdec_SOURCES=\
	Utilities.cpp \
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

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
