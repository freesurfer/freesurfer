## 
# Makefile.am
##

if USE_LOCAL_TIFF
TIFFCFLAGS=-I$(top_srcdir)/tiff
LIBTIFF=$(top_builddir)/tiff/libtiff.a
else
TIFFCFLAGS=$(TIFF_CFLAGS)
LIBTIFF=$(LIB_TIFF)
endif

AM_CPPFLAGS=\
	-I$(top_srcdir)/glut \
	-I$(top_srcdir)/include \
	-I$(top_srcdir)/dicom \
	-I$(top_srcdir)/unix \
	$(GL_CFLAGS) $(TIFFCFLAGS)

AM_LDFLAGS=$(X_LIBS) $(TIFF_LIBS)

if HAVE_OPENGL_LIBS
bin_PROGRAMS = mris2rgb
mris2rgb_SOURCES=mris2rgb.c
mris2rgb_LDADD= $(addprefix $(top_builddir)/,$(LIBS_MGH)) \
	$(LIBS_OPENGL) -lGL $(LIBTIFF) \
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
