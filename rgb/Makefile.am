#
# Makefile.am for libutils.a
#
# Makefile.am specify i addtion to the standard
AM_CPPFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=

lib_LIBRARIES= librgb.a
librgb_a_SOURCES = $(wildcard *.c)
am_librgb_a_OBJECTS = $(librgb_a_SOURCES:%.c=%.o)

.c.o:
	source='$<' object='$@' libtool=no \
	depfile='$(DEPDIR)/$*.Po' tmpdepfile='$(DEPDIR)/$*.TPo' \
	$(CCDEPMODE) $(depcomp) \
	$(COMPILE) -c `test -f '$<' || echo '$(srcdir)/'`$<

