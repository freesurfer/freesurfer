## 
# Makefile.am
##

if ENABLE_XAWPLUS_APPS
AM_CPPFLAGS=$(GLUT_CFLAGS) \
	-DHAVE_XAWPLUS $(XAWPLUS_CFLAGS) \
	-I$(top_srcdir)/include \
	-I$(top_srcdir)/unix 
else
AM_CPPFLAGS=$(GLUT_CFLAGS) \
	-I$(top_srcdir)/include \
	-I$(top_srcdir)/unix 
endif
AM_LDFLAGS= 

if HAVE_GLUT_LIBS
bin_PROGRAMS = nmovie
nmovie_SOURCES=nmovie.c
if ENABLE_XAWPLUS_APPS
nmovie_LDADD= $(addprefix $(top_srcdir)/,$(LIBS_MGH)) \
	$(GLUT_LIBS) -lXt -lXext -lXaw $(XAWPLUS_LIBS) -lXawPlus
else
nmovie_LDADD= $(addprefix $(top_srcdir)/,$(LIBS_MGH)) \
	$(GLUT_LIBS) -lXt -lXext -lXaw
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
