## 
# Makefile.am - specify in addition to the standard
##

AM_CPPFLAGS=$(GLUT_CFLAGS) $(TCL_CFLAGS) $(tcl_includes) \
	-I$(top_srcdir)/include \
	-I$(top_srcdir)/dicom \
	-fno-strict-aliasing
AM_LDFLAGS=$(LDFLAGS_TCL)

bin_PROGRAMS = tkmedit
tkmedit_SOURCES=tkmedit.c \
	tkmDisplayArea.c \
	tkmFunctionalVolume.c \
	tkmMeditWindow.c \
	tkmedit.h \
	tkmDisplayArea.h \
	tkmFunctionalVolume.h \
	tkmMeditWindow.h \
	 $(top_srcdir)/fsgdf/fsgdf_wrap.c
tkmedit_LDADD= $(addprefix $(top_srcdir)/,$(LIBS_MGH)) \
	$(top_srcdir)/fsgdf/libfsgdf.a \
	$(LIBS_TCL_OPENGL) \
	$(X_LIBS) $(GLUT_LIBS)
tkmedit_LDFLAGS=


# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
