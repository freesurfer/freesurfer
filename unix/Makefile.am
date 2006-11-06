#
# Makefile.am
#

AM_CPPFLAGS=-I$(top_srcdir)/include $(GL_CFLAGS) 
AM_LDFLAGS=

noinst_LIBRARIES= libunix.a
libunix_a_SOURCES = $(wildcard *.c) $(wildcard *.h)
am_libunix_a_OBJECTS = $(libunix_a_SOURCES:%.c=%.o)

#lib_LTLIBRARIES= libunix.la
#libunix_la_SOURCES = $(wildcard *.c) $(wildcard *.h)

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
