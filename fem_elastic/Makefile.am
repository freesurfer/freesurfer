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
AM_LDFLAGS=$(OS_LDFLAGS) -framework Accelerate
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

bin_PROGRAMS=createMorph applyMorph exportGcam surf2vol

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

surf2vol_SOURCES=\
	$(common3d_SOURCES) \
	fsurf2vol.cpp cstats.h surf_types.h \
	surf_powell.cpp surf_powell.h \
	surf_energy.cpp surf_energy.h \
	misc_maths.cpp \
	transformUtils.cpp \
	surf_utils.cpp \
	ZlibStringCompressor.cpp \
	pbCluster_mesh_crop.cpp pbCluster_mesh_crop.h mesh_bfs.h \
	untangler.cpp untangler.h pbmesh_crop.h

if HAVE_MAC_OSX
surf2vol_LDADD=\
	$(addprefix $(top_builddir)/, $(LIBS_MGH)) \
	$(top_builddir)/tetgen/libtetgen.a \
  $(PETSC_LIBS) -lpetscts -lpetscsnes -lpetscksp -lpetscdm -lpetscmat \
  -lpetscvec -lpetsc -lmpich -lpmpich
else
surf2vol_LDADD=\
    $(addprefix $(top_builddir)/, $(LIBS_MGH)) \
    $(top_builddir)/tetgen/libtetgen.a \
    $(PETSC_LIBS) -lpetscts -lpetscsnes -lpetscksp -lpetscdm -lpetscmat \
    -lpetscvec -lpetsc -lmpich \
		$(LIB_GFORTRAN) $(LIB_LAPACK_A) $(LIB_BLAS_A) $(LIB_G2C_A) $(LIB_GFORTRAN)
endif


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
