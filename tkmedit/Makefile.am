## 
# Makefile.am
##

AM_CPPFLAGS=$(GLUT_CFLAGS) $(TCL_CFLAGS) $(tcl_includes) \
	-I$(top_srcdir)/include \
	-fno-strict-aliasing
AM_LDFLAGS=$(LDFLAGS_TCL)

if ENABLE_TCLTK_APPS
bin_PROGRAMS=tkmedit
tkmedit_SOURCES=tkmedit.c \
	tkmDisplayArea.c \
	tkmFunctionalVolume.c \
	tkmMeditWindow.c \
	tkmedit.h \
	tkmDisplayArea.h \
	tkmFunctionalVolume.h \
	tkmMeditWindow.h \
	$(top_srcdir)/fsgdf/fsgdf_wrap.c
tkmedit_LDADD= 	$(top_srcdir)/fsgdf/libfsgdf.a \
	$(addprefix $(top_srcdir)/,$(LIBS_MGH)) \
	$(LIBS_TCL_OPENGL) \
	$(X_LIBS) $(GLUT_LIBS)
tkmedit_LDFLAGS=
else
bin_PROGRAMS=
endif

if ENABLE_TCLTK_APPS
# put a wrapper around the bin, used to setup tcltktixblt enviro vars
install-exec-hook:
	mv $(srcdir)/tkmedit $(bindir)/tkmedit.bin
	echo "#!/bin/tcsh -ef" > $(bindir)/tkmedit
	echo "source \$$FREESURFER_HOME/bin/tcl_setup" >> $(bindir)/tkmedit
	echo "tkmedit.bin \$$argv" >> $(bindir)/tkmedit
	chmod a+x $(bindir)/tkmedit
else
install-exec-hook:
endif

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
