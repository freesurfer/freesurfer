#
# Makefile.am for liblineprof.a
#

#if USE_LOCAL_TIFF
#TIFF_THINGS=-I$(top_srcdir)/tiff -L$(top_srcdir)/tiff
#else
#TIFF_THINGS=$(TIFF_CFLAGS)
#endif
#
#if USE_LOCAL_JPEG
#JPEG_THINGS=-I$(top_srcdir)/jpeg -L$(top_srcdir)/jpeg
#else
#JPEG_THINGS=$(JPEG_CFLAGS)
#endif
#
#if USE_LOCAL_EXPAT
#EXPAT_THINGS=-I$(top_srcdir)/expat -L$(top_srcdir)/expat
#else
#EXPAT_THINGS=$(EXPAT_CFLAGS)
#endif
#
#	$(TIFF_THINGS) $(JPEG_THINGS) $(EXPAT_THINGS) \
#	$(GL_CFLAGS) \
#	$(TCL_CFLAGS) \
#	$(VXL_CFLAGS) \
#	-I$(top_srcdir)/include 

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
TESTS=testLineProf

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
  $(LIB_LAPACK_A) $(LIB_BLAS_A) -lgfortran

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
