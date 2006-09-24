## 
# Makefile.am - specify in addition to the standard
##

AM_CPPFLAGS=$(GLUT_CFLAGS) $(TCL_CFLAGS) $(tcl_includes) \
	-I$(top_srcdir)/include \
	-I$(top_srcdir)/dicom \
	-fno-strict-aliasing
AM_LDFLAGS=$(LDFLAGS_TCL)

if ENABLE_TCLTK_APPS
bin_PROGRAMS=tkmedit
else
bin_PROGRAMS=
endif
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

# put a wrapper around the bin, used to setup tcltktixblt enviro vars
install-exec-hook:
	mv $(bindir)/tkmedit $(bindir)/tkmedit.bin
	cp $(top_srcdir)/scripts/tcltktixblt_wrapper $(bindir)/tkmedit
	echo "tkmedit.bin \$$argv" >> $(bindir)/tkmedit
	chmod a+x $(bindir)/tkmedit

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
