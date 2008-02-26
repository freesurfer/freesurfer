## 
# Makefile.am
##

AM_CPPFLAGS=-I$(top_srcdir)/glut \
	-I$(top_srcdir)/include \
	-I$(top_srcdir)/dicom \
	-I$(top_srcdir)/unix \
	$(X_CFLAGS)

AM_LDFLAGS= \
	$(X_LIBS)

if HAVE_OPENGL_LIBS
bin_PROGRAMS = mris2rgb
mris2rgb_SOURCES=mris2rgb.c
mris2rgb_LDADD= $(addprefix $(top_builddir)/,$(LIBS_MGH)) \
	$(LIBS_OPENGL) -lGL \
	$(top_builddir)/glut/libglut.a
mris2rgb_LDFLAGS=
else
bin_PROGRAMS=
endif

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
