#
# Makefile.am
#
# Makefile.am specify i addtion to the standard
AM_CPPFLAGS=$(X_CFLAGS) $(GLUT_CFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=

noinst_LIBRARIES= liblog.a
liblog_a_SOURCES = $(wildcard *.c)
am_liblog_a_OBJECTS = $(liblog_a_SOURCES:%.c=%.o)

.c.o:
	source='$<' object='$@' libtool=no \
	depfile='$(DEPDIR)/$*.Po' tmpdepfile='$(DEPDIR)/$*.TPo' \
	$(CCDEPMODE) $(depcomp) \
	$(COMPILE) -c `test -f '$<' || echo '$(srcdir)/'`$<



# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
