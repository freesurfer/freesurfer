#
# Makefile.am for libutils.a
#

NrrdIO_dir = NrrdIO
cephes_dir = cephes

if ENABLE_ITK_APPS
ITK_THINGS=$(NIFTI_CFLAGS) -DHAVE_ITK_LIBS
else
ITK_THINGS=
endif

if HAVE_NIFTI_LIBS
NIFTI_THINGS=$(NIFTI_CFLAGS) -DHAVE_NIFTI_LIBS
else
NIFTI_THINGS=
endif

if HAVE_GLUT_LIBS
GLUT_THINGS=$(GLUT_CFLAGS) -DHAVE_GLUT_LIBS
else
GLUT_THINGS=
endif

AM_CPPFLAGS=$(X_CFLAGS) $(VXL_CFLAGS) \
	-I$(top_srcdir)/include \
	-I$(top_srcdir)/include/dicom \
	-I$(top_srcdir)/include/NrrdIO \
	$(ITK_THINGS) $(NIFTI_THINGS) $(GLUT_THINGS) -DHAVE_ZLIB

AM_LDFLAGS=

noinst_LIBRARIES= libutils.a

libutils_a_SOURCES = $(wildcard *.c) $(wildcard *.cpp)

am_libutils_a_OBJECTS = $(libutils_a_SOURCES:%.c=%.o) \
	$(libutils_a_SOURCES:%.cpp=%.o) \
	$(wildcard $(cephes_dir)/*.o) \
	$(wildcard $(NrrdIO_dir)/*.o)

#lib_LTLIBRARIES= libutils.la
#libutils_la_SOURCES = $(wildcard *.c) $(wildcard *.h)

if HAVE_MAC_OSX
# mrisurf must be compiled w/o optimization enabled, in order 
# to mask a segfault appearing on the Mac when running mris_fix_topology
mrisurf.o: mrisurf.c
	@echo "compiling mrisurf.c to mrisurf.o (with -O0)"
	source='$<' object='$@' libtool=no \
	depfile='$(DEPDIR)/$*.Po' tmpdepfile='$(DEPDIR)/$*.TPo' \
	$(CCDEPMODE) $(depcomp) \
	$(COMPILE) -O0 -c `test -f '$<' || echo '$(srcdir)/'`$<
endif

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
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
