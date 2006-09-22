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

bin_PROGRAMS = tksurfer
tksurfer_SOURCES=tksurfer.c $(top_srcdir)/fsgdf/fsgdf_wrap.c
tksurfer_LDADD= $(addprefix $(top_srcdir)/,$(LIBS_MGH)) \
	$(top_srcdir)/fsgdf/libfsgdf.a \
	$(LIBS_TCL_OPENGL) \
	$(X_LIBS) $(GLUT_LIBS) -lGL
tksurfer_LDFLAGS=

install-exec-hook:
        mv $(bindir)/tksurfer $(bindir)/tksurfer.bin
        cp tksurfer_wrapper $(bindir)/tksurfer
        chmod a+x $(bindir)/tksurfer

EXTRA_DIST="tksurfer_wrapper"

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
