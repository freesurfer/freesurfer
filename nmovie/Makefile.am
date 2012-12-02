## 
# Makefile.am
##

if ENABLE_XAWPLUS_APPS
AM_CPPFLAGS=\
	-DHAVE_XAWPLUS $(XAWPLUS_CFLAGS) \
	-DHAVE_OPENGL $(GL_CFLAGS) \
	-I$(top_srcdir)/include \
	-I$(top_srcdir)/unix \
	-I$(top_srcdir)/glut
else
AM_CPPFLAGS=\
	-DHAVE_OPENGL $(GL_CFLAGS) \
	-I$(top_srcdir)/include \
	-I$(top_srcdir)/unix \
	-I$(top_srcdir)/glut
endif
AM_LDFLAGS= 

if HAVE_OPENGL_LIBS
bin_PROGRAMS = nmovie
nmovie_SOURCES=nmovie.c
if ENABLE_XAWPLUS_APPS
nmovie_LDADD= $(addprefix $(top_builddir)/,$(LIBS_MGH)) \
	$(LIBS_OPENGL) $(XAWPLUS_LIBS) $(LIBS_XAWPLUS) -lXext \
	$(top_builddir)/glut/libglut.a
else
nmovie_LDADD= $(addprefix $(top_builddir)/,$(LIBS_MGH)) \
	$(LIBS_OPENGL) -lXt -lXext -lXaw \
	$(top_builddir)/glut/libglut.a
endif
nmovie_LDFLAGS=$(X_LIBS)
else
bin_PROGRAMS=
endif

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
