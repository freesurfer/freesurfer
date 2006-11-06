## 
# Makefile.am
##

AM_CPPFLAGS=$(GLUT_CFLAGS) \
	-I$(top_srcdir)/include \
	-I$(top_srcdir)/unix
AM_LDFLAGS=

if HAVE_GLUT_LIBS
bin_PROGRAMS = mris_show
mris_show_SOURCES=mris_show.c
mris_show_LDADD= $(addprefix $(top_srcdir)/,$(LIBS_MGH)) $(GLUT_LIBS)
mris_show_LDFLAGS=$(X_LIBS)
else
bin_PROGRAMS=
endif

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
