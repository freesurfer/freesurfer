#
# Makefile.am for liblineprof.a
#

# For now skip the lineprof build when --disable-all-apps is set (hence, 
# ENABLE_QT_APPS = false) because lineprof is only used by freeview.
if ENABLE_QT_APPS

if USE_LOCAL_TIFF
TIFFCFLAGS=-I$(top_srcdir)/tiff
LIBTIFF=$(top_builddir)/tiff/libtiff.a
else
TIFFCFLAGS=$(TIFF_CFLAGS)
LIBTIFF=$(LIB_TIFF)
endif

AM_CPPFLAGS= -Wno-deprecated \
	-I$(top_srcdir)/include \
	-I$(top_srcdir)/tetgen \
	$(VXL_CFLAGS) \
	$(VTK_CXXFLAGS) \
  $(ITK_CFLAGS) -DHAVE_ITK_LIBS \
  $(PETSC_CFLAGS) $(TIFFCFLAGS) \
	$(BOOST_CFLAGS)

if HAVE_MAC_OSX
AM_LDFLAGS=$(OS_LDFLAGS) -framework Accelerate -framework Cocoa
else
AM_LDFLAGS=$(OS_LDFLAGS) -z muldefs
endif

AM_CXXFLAGS=$(AM_CPPFLAGS) $(VTK_COCOA_CXXFLAGS)

noinst_LIBRARIES= liblineprof.a

liblineprof_a_SOURCES= \
	Tracer.cpp \
	Preprocessor.cpp \
	PetscSolver.cpp \
	LineProf.cpp

#foo_SOURCES=
#foo:
#	if ! test -f "$(top_builddir)/vtkutils/vtkKWRGBATransferFunctionEditorTester-scalars.mgh"; then \
#	cp -v $(top_srcdir)/vtkutils/vtkKWRGBATransferFunctionEditorTester-scalars.mgh \
#				$(top_builddir)/vtkutils/; fi;
#	if ! test -f "$(top_builddir)/vtkutils/IconLoaderTestIcons.txt"; then \
#	cp -v $(top_srcdir)/vtkutils/IconLoaderTestIcons.txt \
#				$(top_builddir)/vtkutils/; fi;
##check_PROGRAMS=foo LineProfTest

check_PROGRAMS=LineProfTest
#NOTE: the testLineProf script calls the vtk setup script, which is
# necessary under normal circumstances, but since configure sets the
# TEST_ENVIRONMENT var already, it is not needed for make check
# (and in fact will not work for make distcheck), so LineProfTest
# can be called directly
#TESTS=testLineProf
TESTS=LineProfTest

LineProfTest_SOURCES=LineProfTest.cpp


if HAVE_MAC_OSX
LineProfTest_LDADD=\
	$(top_builddir)/lineprof/liblineprof.a \
	$(VTK_COCOA_LIBS) $(VTK_VERDICT_LIB) $(VTK_SQLITE_LIB) \
	-lvtkIO -lvtkWidgets -lvtkGraphics \
	-lvtkImaging -lvtkFiltering \
	-lvtkCommon -lvtksys \
	-lvtkHybrid -lvtkexoIIc -lvtkNetCDF -lvtkDICOMParser \
	-lvtkVolumeRendering -lvtkRendering -lvtkftgl \
	-ltiff $(VTK_TIFF_LIB) -ljpeg $(VTK_JPEG_LIB) \
	-lexpat $(VTK_EXPAT_LIB) $(VTK_FREETYPE_LIB) \
	$(VTK_PNG_LIB) $(VTK_Z_LIB) $(VTK_METAIO_LIB) \
  $(ITK_LIBS) \
  $(PETSC_LIBS) -lpetscts -lpetscsnes -lpetscksp \
  -lpetscdm -lpetscmat -lpetscvec -lpetsc -lmpich -lpmpich
else
LineProfTest_LDADD=\
	$(top_builddir)/lineprof/liblineprof.a \
	$(VTK_LIBS) $(VTK_VERDICT_LIB) $(VTK_SQLITE_LIB) \
	-lvtkIO -lvtkWidgets -lvtkGraphics \
	-lvtkImaging -lvtkFiltering \
	-lvtkCommon -lvtksys \
	-lvtkHybrid -lvtkexoIIc -lvtkNetCDF -lvtkDICOMParser \
	-lvtkVolumeRendering -lvtkRendering -lvtkftgl \
	-lexpat $(VTK_EXPAT_LIB) $(VTK_FREETYPE_LIB) \
    $(LIBTIFF) $(LIBJPEG) \
	$(VTK_PNG_LIB) $(VTK_Z_LIB) $(VTK_METAIO_LIB) \
  $(ITK_LIBS) \
  $(PETSC_LIBS) -lpetscts -lpetscsnes -lpetscksp \
  -lpetscdm -lpetscmat -lpetscvec -lpetsc -lmpich -lfmpich \
  $(LIB_GFORTRAN) $(LIB_LAPACK_A) $(LIB_BLAS_A) $(LIB_G2C_A) $(LIB_GFORTRAN)
endif

endif
# endif ENABLE_QT_APPS

EXTRA_DIST=  \
PetscSolver.h \
Preprocessor.h \
Tracer.h \
LineProf.h \
LineProfTest.cpp \
testLineProf

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
