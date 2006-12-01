#
# Makefile.am
#
AM_CPPFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom -fno-strict-aliasing
AM_LDFLAGS=

noinst_LIBRARIES= libdicom.a
libdicom_a_SOURCES = $(wildcard *.c)
am_libdicom_a_OBJECTS = $(libdicom_a_SOURCES:%.c=%.o)

.c.o:
	source='$<' object='$@' libtool=no \
	depfile='$(DEPDIR)/$*.Po' tmpdepfile='$(DEPDIR)/$*.TPo' \
	$(CCDEPMODE) $(depcomp) \
	$(COMPILE) -c `test -f '$<' || echo '$(srcdir)/'`$<

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
