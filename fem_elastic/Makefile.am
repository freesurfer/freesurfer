#
# Makefile.am
#

if HAVE_ITK_LIBS
if HAVE_PETSC_LIBS
if HAVE_BOOST_LIBS
AM_CPPFLAGS=\
	-I$(top_srcdir)/include \
	-I$(top_srcdir)/tetgen \
	$(VXL_CFLAGS)	\
	$(ITK_CFLAGS) -DHAVE_ITK_LIBS \
	$(PETSC_CFLAGS) \
	$(BOOST_CFLAGS)

if HAVE_MAC_OSX
AM_LDFLAGS=$(OS_LDFLAGS)
else
AM_LDFLAGS=$(OS_LDFLAGS) -z muldefs
endif

common3d_SOURCES=\
	fem_3d.cpp fem_3d.h solver.h \
	misc.h mesh.h toctree.hpp node.h simple_timer.h tag_fio.h \
	coords.h material.h element.h gmpError.h streamIoMisc.h stl_utils.hpp \
	small_matrix.cpp small_matrix.h \
	morph.cpp morph.h \
	morph_utils.cpp morph_utils.h

bin_PROGRAMS=createMorph applyMorph exportGcam

createMorph_SOURCES=\
	$(common3d_SOURCES) \
	fcreateMorph.cpp \
	misc_maths.cpp misc_maths.h \
	transformUtils.cpp transformUtils.h \
	surf_utils.cpp surf_utils.h \
	ZlibStringCompressor.cpp ZlibStringCompressor.h

createMorph_LDADD=\
	$(addprefix $(top_builddir)/, $(LIBS_MGH)) \
	$(top_builddir)/tetgen/libtetgen.a $(BOOST_PO_LIB)

applyMorph_SOURCES=\
	$(common3d_SOURCES) \
	fapplyMorph.cpp \
	misc_maths.cpp misc_maths.h \
	transformUtils.cpp transformUtils.h \
	surf_utils.cpp surf_utils.h \
	ZlibStringCompressor.cpp ZlibStringCompressor.h

applyMorph_LDADD=\
	$(addprefix $(top_builddir)/, $(LIBS_MGH)) \
	$(top_builddir)/tetgen/libtetgen.a $(BOOST_PO_LIB)

exportGcam_SOURCES=\
	$(common3d_SOURCES) \
	exportGcam.cpp \
	ZlibStringCompressor.cpp ZlibStringCompressor.h

exportGcam_LDADD=\
	$(addprefix $(top_builddir)/, $(LIBS_MGH)) \
	$(top_builddir)/tetgen/libtetgen.a $(BOOST_PO_LIB)

# build library...
noinst_LIBRARIES= libfem_elastic.a
libfem_elastic_a_SOURCES = $(createMorph_SOURCES) $(applyMorph_SOURCES)

endif
endif
endif

EXTRA_DIST=

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
