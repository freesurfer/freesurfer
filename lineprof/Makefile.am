#
# Makefile.am for liblineprof.a
#

AM_CPPFLAGS= -Wno-deprecated \
	$(VXL_CFLAGS) \
	$(VTK_CXXFLAGS) \
  $(ITK_CFLAGS) \
  $(PETSC_CFLAGS)

AM_CXXFLAGS=$(AM_CPPFLAGS)

AM_LDFLAGS=


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

LineProfTest_LDADD=\
	$(top_builddir)/lineprof/liblineprof.a \
	$(VTK_LIBS) $(VTK_VERDICT_LIB) $(VTK_SQLITE_LIB) \
	-lvtkIO -lvtkWidgets -lvtkGraphics \
	-lvtkImaging -lvtkFiltering \
	-lvtkCommon -lvtksys \
	-lvtkHybrid -lvtkexoIIc -lvtkNetCDF -lvtkDICOMParser \
	-lvtkVolumeRendering -lvtkRendering -lvtkftgl \
	-lexpat $(VTK_EXPAT_LIB) $(VTK_FREETYPE_LIB) \
	-ltiff $(VTK_TIFF_LIB) -ljpeg $(VTK_JPEG_LIB) \
	$(VTK_PNG_LIB) $(VTK_Z_LIB) $(VTK_METAIO_LIB) \
  $(ITK_LIBS) \
  $(PETSC_LIBS) -lpetscts -lpetscsnes -lpetscksp \
  -lpetscdm -lpetscmat -lpetscvec -lpetsc -lmpich -lfmpich \
  $(LIB_GFORTRAN) $(LIB_LAPACK_A) $(LIB_BLAS_A) $(LIB_G2C_A) $(LIB_GFORTRAN)

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
