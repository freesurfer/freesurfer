#
# Makefile.am
#

AM_CPPFLAGS=-I$(top_srcdir)/include \
	$(GL_CFLAGS) \
	$(HAVE_OPENGL)

AM_LDFLAGS=

noinst_LIBRARIES=libunix.a

libunix_a_SOURCES=cuserid.c unix_thread.c xwin.c \
	ftime.c windiag.c xwindow.c \
	oglutil.c xvmri.c texfont.c xvutil.c glStub.c

EXTRA_DIST=lwpproto.h xvmri.h xwindow.h \
	oglutil.h xvutil.h xwin.h

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
