## 
# Makefile.am
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
	tkmMeditWindow.h
tkmedit_LDADD= $(addprefix $(top_srcdir)/,$(LIBS_MGH)) \
	$(LIBS_TCL_OPENGL) \
	$(X_LIBS) $(GLUT_LIBS)
tkmedit_LDFLAGS=

install-exec-hook:
        mv $(bindir)/tkmedit $(bindir)/tkmedit.bin
        cp tkmedit_wrapper $(bindir)/tkmedit
        chmod a+x $(bindir)/tkmedit

EXTRA_DIST="tkmedit_wrapper"

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
