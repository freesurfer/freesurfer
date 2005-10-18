#
# Makefile.am for libutils.a
#

AM_CPPFLAGS=$(X_CFLAGS) -I$(top_srcdir)/include -I$(top_srcdir)/include/dicom $(GSL_CFLAGS)
AM_LDFLAGS=$(GSL_LIBS)

#lib_LTLIBRARIES = libutils.la
#libutils_la_SOURCES = $(wildcard *.c)
#libutils_la_LDFLAGS = version-info 0:0:0

noinst_LIBRARIES= libutils.a
libutils_a_SOURCES = $(wildcard *.c) $(wildcard *.h)
am_libutils_a_OBJECTS = $(libutils_a_SOURCES:%.c=%.o)

# tqli needs with -ffloat-store
tqli.o: tqli.c
	@echo "compiling tqli.c to tqli.o (with -ffloat-store)"
	source='$<' object='$@' libtool=no \
	depfile='$(DEPDIR)/$*.Po' tmpdepfile='$(DEPDIR)/$*.TPo' \
	$(CCDEPMODE) $(depcomp) \
	$(COMPILE) -ffloat-store -c `test -f '$<' || echo '$(srcdir)/'`$<

# mrisurf must be compiled w/o optimization disabled, in order 
# to mask a segfault appearing on the Mac when running mris_fix_topology
mrisurf.o: mrisurf.c
	@echo "compiling mrisurf.c to mrisurf.o (with -O0)"
	source='$<' object='$@' libtool=no \
	depfile='$(DEPDIR)/$*.Po' tmpdepfile='$(DEPDIR)/$*.TPo' \
	$(CCDEPMODE) $(depcomp) \
	$(COMPILE) -O0 -c `test -f '$<' || echo '$(srcdir)/'`$<

.c.o:
	source='$<' object='$@' libtool=no \
	depfile='$(DEPDIR)/$*.Po' tmpdepfile='$(DEPDIR)/$*.TPo' \
	$(CCDEPMODE) $(depcomp) \
	$(COMPILE) -c `test -f '$<' || echo '$(srcdir)/'`$<

SUBDIRS=test
