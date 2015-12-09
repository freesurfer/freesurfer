## 
# Makefile.am
##

AM_CPPFLAGS=\
	-I$(top_srcdir)/glut \
	-I$(top_srcdir)/include \
	-I$(top_srcdir)/unix $(GL_CFLAGS)
AM_LDFLAGS=

if HAVE_OPENGL_LIBS
bin_PROGRAMS = mris_show
mris_show_SOURCES=mris_show.c
mris_show_LDADD= $(addprefix $(top_builddir)/,$(LIBS_MGH)) \
	$(LIBS_OPENGL) $(top_builddir)/glut/libglut.a $(X_LIBS)
mris_show_LDFLAGS=
else
bin_PROGRAMS=
endif

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
