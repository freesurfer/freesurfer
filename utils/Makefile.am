#
# Makefile.am for libutils.a
#

NrrdIO_dir = NrrdIO

AM_CPPFLAGS=$(X_CFLAGS) \
	-I$(top_srcdir)/include \
	-I$(top_srcdir)/include/dicom \
	-I$(top_srcdir)/include/NrrdIO \
	$(GSL_CFLAGS) -DHAVE_ZLIB

AM_LDFLAGS=$(GSL_LIBS)

#lib_LTLIBRARIES = libutils.la
#libutils_la_SOURCES = $(wildcard *.c)
#libutils_la_LDFLAGS = version-info 0:0:0

noinst_LIBRARIES= libutils.a
libutils_a_SOURCES = $(wildcard *.c)
am_libutils_a_OBJECTS = $(libutils_a_SOURCES:%.c=%.o) \
	754.o mop.o array.o parseAir.o dio.o sane.o endianAir.o string.o enum.o miscAir.o biffbiff.o accessors.o defaultsNrrd.o enumsNrrd.o arraysNrrd.o methodsNrrd.o reorder.o axis.o simple.o comment.o keyvalue.o endianNrrd.o parseNrrd.o gzio.o read.o write.o format.o formatNRRD.o encoding.o encodingRaw.o encodingAscii.o encodingHex.o encodingGzip.o subset.o encodingBzip2.o formatEPS.o formatPNG.o formatPNM.o formatText.o formatVTK.o


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

SUBDIRS=$(NrrdIO_dir) test

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
