## 
# Makefile.am -specify in addition to the standard
##

AM_CPPFLAGS=$(GLUT_CFLAGS) \
	$(X_CFLAGS) \
	-I$(top_srcdir)/include \
	-I$(top_srcdir)/dicom \
	-I$(top_srcdir)/unix
AM_LDFLAGS= \
	$(X_LIBS)

bin_PROGRAMS = mris2rgb
mris2rgb_SOURCES=mris2rgb.c
mris2rgb_LDADD= $(addprefix $(top_srcdir)/,$(LIBS_MGH)) \
	$(GLUT_LIBS) -lGL
mris2rgb_LDFLAGS=


# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
