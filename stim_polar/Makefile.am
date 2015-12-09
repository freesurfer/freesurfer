## 
# Makefile.am
##

AM_CFLAGS=\
	-I$(top_srcdir)/include \
	-I$(top_srcdir)/glut $(GL_CFLAGS)
AM_LDFLAGS=

if HAVE_OPENGL_LIBS
bin_PROGRAMS = stim_polar
stim_polar_SOURCES=stim_polar.c $(top_srcdir)/utils/machine.c
stim_polar_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH)) \
	$(LIBS_OPENGL) $(top_builddir)/glut/libglut.a $(X_LIBS) 
stim_polar_LDFLAGS=
else
bin_PROGRAMS=
endif

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
