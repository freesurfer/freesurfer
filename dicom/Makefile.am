#
# Makefile.am
#
# Makefile.am specify i addtion to the standard
AM_CPPFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=

noinst_LIBRARIES= libdicom.a
libdicom_a_SOURCES = $(wildcard *.c)
am_libdicom_a_OBJECTS = $(libdicom_a_SOURCES:%.c=%.o)

.c.o:
	source='$<' object='$@' libtool=no \
	depfile='$(DEPDIR)/$*.Po' tmpdepfile='$(DEPDIR)/$*.TPo' \
	$(CCDEPMODE) $(depcomp) \
	$(COMPILE) -c `test -f '$<' || echo '$(srcdir)/'`$<

