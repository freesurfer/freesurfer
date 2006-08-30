#
# Makefile.am for libutils.a
#

NrrdIO_dir = NrrdIO
cephes_dir = cephes

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

libutils_a_SOURCES = $(wildcard *.c) $(wildcard *.cpp)

am_libutils_a_OBJECTS = $(libutils_a_SOURCES:%.c=%.o) \
	$(libutils_a_SOURCES:%.cpp=%.o) \
	$(NrrdIO_dir)/754.o \
	$(NrrdIO_dir)/mop.o \
	$(NrrdIO_dir)/array.o \
	$(NrrdIO_dir)/parseAir.o \
	$(NrrdIO_dir)/dio.o \
	$(NrrdIO_dir)/sane.o \
	$(NrrdIO_dir)/endianAir.o \
	$(NrrdIO_dir)/string.o \
	$(NrrdIO_dir)/enum.o \
	$(NrrdIO_dir)/miscAir.o \
	$(NrrdIO_dir)/biffbiff.o \
	$(NrrdIO_dir)/accessors.o \
	$(NrrdIO_dir)/defaultsNrrd.o \
	$(NrrdIO_dir)/enumsNrrd.o \
	$(NrrdIO_dir)/arraysNrrd.o \
	$(NrrdIO_dir)/methodsNrrd.o \
	$(NrrdIO_dir)/reorder.o \
	$(NrrdIO_dir)/axis.o \
	$(NrrdIO_dir)/simple.o \
	$(NrrdIO_dir)/comment.o \
	$(NrrdIO_dir)/keyvalue.o \
	$(NrrdIO_dir)/endianNrrd.o \
	$(NrrdIO_dir)/parseNrrd.o \
	$(NrrdIO_dir)/gzio.o \
	$(NrrdIO_dir)/read.o \
	$(NrrdIO_dir)/write.o \
	$(NrrdIO_dir)/format.o \
	$(NrrdIO_dir)/formatNRRD.o \
	$(NrrdIO_dir)/encoding.o \
	$(NrrdIO_dir)/encodingRaw.o \
	$(NrrdIO_dir)/encodingAscii.o \
	$(NrrdIO_dir)/encodingHex.o \
	$(NrrdIO_dir)/encodingGzip.o \
	$(NrrdIO_dir)/subset.o \
	$(NrrdIO_dir)/encodingBzip2.o \
	$(NrrdIO_dir)/formatEPS.o \
	$(NrrdIO_dir)/formatPNG.o \
	$(NrrdIO_dir)/formatPNM.o \
	$(NrrdIO_dir)/formatText.o \
	$(NrrdIO_dir)/formatVTK.o \
	$(cephes_dir)/bdtr.o \
	$(cephes_dir)/chdtr.o \
	$(cephes_dir)/const.o \
	$(cephes_dir)/fdtr.o \
	$(cephes_dir)/gamma.o \
	$(cephes_dir)/igam.o \
	$(cephes_dir)/igami.o \
	$(cephes_dir)/incbet.o \
	$(cephes_dir)/incbi.o \
	$(cephes_dir)/mtherr.o \
	$(cephes_dir)/ndtr.o \
	$(cephes_dir)/ndtri.o \
	$(cephes_dir)/polevl.o \
	$(cephes_dir)/stdtr.o 

# mrisurf must be compiled w/o optimization enabled, in order 
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

.cpp.o:
	@echo "compiler: $(COMPILE) end_compiler"
	source='$<' object='$@' libtool=no \
	depfile='$(DEPDIR)/$*.Po' tmpdepfile='$(DEPDIR)/$*.TPo' \
	$(CCDEPMODE) $(depcomp) \
	$(COMPILE) -c `test -f '$<' || echo '$(srcdir)/'`$<

SUBDIRS=$(NrrdIO_dir) $(cephes_dir) test

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
