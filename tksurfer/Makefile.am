## 
# Makefile.am
##

AM_CPPFLAGS=$(GLUT_CFLAGS) \
	$(TCL_CFLAGS) \
	$(X_CFLAGS) \
	-I$(top_srcdir)/include \
	-I$(top_srcdir)/dicom \
	-I$(top_srcdir)/unix

AM_LDFLAGS=$(LDFLAGS_TCL)

if ENABLE_TCLTK_APPS
bin_PROGRAMS=tksurfer
tksurfer_SOURCES=tksurfer.c $(top_srcdir)/fsgdf/fsgdf_wrap.c
tksurfer_LDADD= 	$(top_srcdir)/fsgdf/libfsgdf.a \
  $(addprefix $(top_srcdir)/,$(LIBS_MGH)) \
	$(LIBS_TCL_OPENGL) \
	$(X_LIBS) $(GLUT_LIBS) -lGL
tksurfer_LDFLAGS=
else
bin_PROGRAMS=
endif

if ENABLE_TCLTK_APPS
install-exec-hook:
	mv $(srcdir)/tksurfer $(bindir)/tksurfer.bin
	echo "#!/bin/tcsh -ef" > $(bindir)/tksurfer
	echo "source \$$FREESURFER_HOME/bin/tcl_setup" >> $(bindir)/tksurfer
	echo "tksurfer.bin \$$argv" >> $(bindir)/tksurfer
	chmod a+x $(bindir)/tksurfer
else
install-exec-hook:
endif

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
