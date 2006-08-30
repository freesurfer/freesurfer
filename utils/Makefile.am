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
  $(wildcard $(cephes_dir)/*.o) \
	$(wildcard $(NrrdIO_dir)/*.o)

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
