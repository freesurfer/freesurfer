#
# Makefile.am for libqdecproject.a
#

AM_CPPFLAGS= -Wno-deprecated \
	$(X_CFLAGS) \
	$(VXL_CFLAGS) \
	-I$(top_srcdir)/include

AM_LDFLAGS=

noinst_LIBRARIES= libqdecproject.a

libqdecproject_a_SOURCES= QdecContrast.cxx  QdecDataTable.cxx  \
	QdecFactor.cxx QdecGlmDesign.cxx  QdecGlmFit.cxx  QdecGlmFitResults.cxx  \
	QdecProject.cxx  QdecSubject.cxx

am_libqdecprojet_a_OBJECTS=$(libqdecproject_a_SOURCES:%.cxx=%.o)

EXTRA_DIST=	ProgressUpdateGUI.h  QdecContrast.h  QdecDataTable.h  \
	QdecFactor.h  QdecGlmDesign.h  QdecGlmFit.h  QdecGlmFitResults.h  \
	QdecProject.h  QdecSubject.h  umbrello.xmi

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
