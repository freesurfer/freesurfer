#
# Makefile.am for libutils.a
#

AM_CPPFLAGS=$(X_CFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=

#lib_LTLIBRARIES = libutils.la
#libutils_la_SOURCES = $(wildcard *.c)
#libutils_la_LDFLAGS = version-info 0:0:0

noinst_LIBRARIES= libutils.a
libutils_a_SOURCES = $(wildcard *.c)
am_libutils_a_OBJECTS = $(libutils_a_SOURCES:%.c=%.o)

# tqli needs with -ffloat-store
tqli.o: tqli.c
	@echo "compiling tqli.c to tqli.o (with -ffloat-store)"
	source='$<' object='$@' libtool=no \
	depfile='$(DEPDIR)/$*.Po' tmpdepfile='$(DEPDIR)/$*.TPo' \
	$(CCDEPMODE) $(depcomp) \
	$(COMPILE) -ffloat-store -c `test -f '$<' || echo '$(srcdir)/'`$<

.c.o:
	source='$<' object='$@' libtool=no \
	depfile='$(DEPDIR)/$*.Po' tmpdepfile='$(DEPDIR)/$*.TPo' \
	$(CCDEPMODE) $(depcomp) \
	$(COMPILE) -c `test -f '$<' || echo '$(srcdir)/'`$<

SUBDIRS=test
