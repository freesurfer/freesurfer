#
# Makefile.am
#
# Makefile.am specify i addtion to the standard
AM_CPPFLAGS=$(X_CFLAGS) $(GLUT_CFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=

lib_LIBRARIES= libunix.a
libunix_a_SOURCES = $(wildcard *.c)
am_libunix_a_OBJECTS = $(libunix_a_SOURCES:%.c=%.o)

.c.o:
	source='$<' object='$@' libtool=no \
	depfile='$(DEPDIR)/$*.Po' tmpdepfile='$(DEPDIR)/$*.TPo' \
	$(CCDEPMODE) $(depcomp) \
	$(COMPILE) -c `test -f '$<' || echo '$(srcdir)/'`$<

