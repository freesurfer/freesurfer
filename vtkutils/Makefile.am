#
# Makefile.am for libvtkutils.a
#

if USE_LOCAL_TIFF
TIFF_THINGS=-I$(top_srcdir)/tiff -L$(top_srcdir)/tiff
else
TIFF_THINGS=$(TIFF_CFLAGS)
endif

if USE_LOCAL_JPEG
JPEG_THINGS=-I$(top_srcdir)/jpeg -L$(top_srcdir)/jpeg
else
JPEG_THINGS=$(JPEG_CFLAGS)
endif

if USE_LOCAL_EXPAT
EXPAT_THINGS=-I$(top_srcdir)/expat -L$(top_srcdir)/expat
else
EXPAT_THINGS=$(EXPAT_CFLAGS)
endif


AM_CPPFLAGS= -Wno-deprecated \
	$(TIFF_THINGS) $(JPEG_THINGS) $(EXPAT_THINGS) \
	$(GL_CFLAGS) \
	$(TCL_CFLAGS) \
	$(VXL_CFLAGS) \
	-I$(top_srcdir)/include \
	$(VTK_CXXFLAGS) \
	$(KWWIDGETS_CXXFLAGS)

AM_CXXFLAGS=$(AM_CPPFLAGS)

AM_LDFLAGS=

if ENABLE_KWWIDGETS_APPS
noinst_LIBRARIES= libvtkutils.a
else
noinst_LIBRARIES=
endif

libvtkutils_a_SOURCES= \
	vtkKWProgressDialog.cxx \
	vtkKWRGBATransferFunctionEditor.cxx \
	vtkKWRGBATransferFunctionEditorTcl.cxx \
	vtkKWBltGraph.cxx vtkKWBltGraphTcl.cxx \
	vtkFreesurferLookupTable.cxx \
	vtkArrowPipeline.cxx \
	vtkRGBATransferFunction.cxx \
	vtkRGBAColorTransferFunction.cxx \
	vtkFDTensorGlyph.cxx \
	vtkODFGlyph.cxx \
	vtkInflatePolyData.cxx \
	IconLoader.cxx

# We use vtkWrapTcl to create the Tcl wrapped code for this object.
vtkKWRGBATransferFunctionEditorTcl.cxx: vtkKWRGBATransferFunctionEditor.h
	${VTK_DIR}/bin/vtkWrapTcl vtkKWRGBATransferFunctionEditor.h \
	${VTK_LIB_VTK}/hints 1 vtkKWRGBATransferFunctionEditorTcl.cxx

vtkKWBltGraphTcl.cxx: vtkKWBltGraph.h
	${VTK_DIR}/bin/vtkWrapTcl vtkKWBltGraph.h \
	${VTK_LIB_VTK}/hints 1 vtkKWBltGraphTcl.cxx

if ENABLE_KWWIDGETS_APPS
foo_SOURCES=
foo:
	if ! test -f "$(top_builddir)/vtkutils/vtkKWRGBATransferFunctionEditorTester-scalars.mgh"; then \
	cp -v $(top_srcdir)/vtkutils/vtkKWRGBATransferFunctionEditorTester-scalars.mgh \
				$(top_builddir)/vtkutils/; fi;
	if ! test -f "$(top_builddir)/vtkutils/IconLoaderTestIcons.txt"; then \
	cp -v $(top_srcdir)/vtkutils/IconLoaderTestIcons.txt \
				$(top_builddir)/vtkutils/; fi;

if HAVE_MAC_OSX
check_PROGRAMS=foo vtkArrowPipelineTester
else
check_PROGRAMS=foo IconLoaderTest \
vtkKWRGBATransferFunctionEditorTester vtkArrowPipelineTester
endif
#these commented-out because the tests are interactive (can't run automatically)
#TESTS=vtkKWRGBATransferFunctionEditorTester vtkArrowPipelineTester
#TESTS=testIconLoader
else
check_PROGRAMS=
TESTS=
endif

IconLoaderTest_SOURCES=IconLoaderTest.cxx IconLoaderTest.h IconLoader.cxx

IconLoaderTest_LDADD=\
	$(top_builddir)/vtkfsio/libvtkfsio.a \
	$(top_builddir)/vtkutils/libvtkutils.a \
	$(LDFLAGS_TCL) $(LIBS_TCL_OPENGL) \
	$(VTK_LIBS) $(VTK_VERDICT_LIB) $(VTK_SQLITE_LIB) \
	-lvtkIO -lvtkWidgets -lvtkGraphics \
	-lvtkImaging -lvtkFiltering \
	-lvtkCommon -lvtksys \
	-lvtkHybrid -lvtkexoIIc -lvtkNetCDF \
	-lvtkVolumeRendering -lvtkRendering -lvtkftgl \
	-lexpat $(VTK_EXPAT_LIB) $(VTK_FREETYPE_LIB) \
	-ltiff $(VTK_TIFF_LIB) -ljpeg $(VTK_JPEG_LIB) \
	$(VTK_PNG_LIB) $(VTK_Z_LIB) $(VTK_METAIO_LIB) \
	-lvtkIO	-lvtkDICOMParser $(VTK_MPEG2ENCODE_LIB) \
	-lvtkImagingTCL -lvtkVolumeRenderingTCL -lvtkRenderingTCL \
	-lvtkFilteringTCL -lvtkWidgetsTCL \
	-lvtkHybridTCL -lvtkGraphicsTCL -lvtkImagingTCL -lvtkIOTCL \
	-lvtkCommonTCL \
	$(KWWIDGETS_LIBS) -lKWWidgets


# vtkKWRGBATransferFunctionEditorTester test program
vtkKWRGBATransferFunctionEditorTester_SOURCES=vtkKWRGBATransferFunctionEditorTester.cxx \
	vtkKWRGBATransferFunctionEditorTesterTcl.cxx \
	RgbatransferfunctioneditortesterLibInit.cxx

vtkKWRGBATransferFunctionEditorTesterTcl.cxx: vtkKWRGBATransferFunctionEditorTester.h
	${VTK_DIR}/bin/vtkWrapTcl vtkKWRGBATransferFunctionEditorTester.h \
	${VTK_LIB_VTK}/hints 1 vtkKWRGBATransferFunctionEditorTesterTcl.cxx

RgbatransferfunctioneditortesterLibInit.cxx: RgbatransferfunctioneditortesterLibInit.data
	${VTK_DIR}/bin/vtkWrapTclInit RgbatransferfunctioneditortesterLibInit.data \
	RgbatransferfunctioneditortesterLibInit.cxx

vtkKWRGBATransferFunctionEditorTester_LDADD=\
	$(top_builddir)/fsgdf/libfsgdf.a \
	$(top_builddir)/vtkfsio/libvtkfsio.a \
	$(top_builddir)/vtkutils/libvtkutils.a \
	$(addprefix $(top_builddir)/, $(LIBS_MGH)) \
	$(LDFLAGS_TCL) $(LIBS_TCL_OPENGL) \
	$(VTK_LIBS) $(VTK_VERDICT_LIB) $(VTK_SQLITE_LIB) \
	-lvtkWidgets -lvtkGraphics \
	-lvtkImaging -lvtkFiltering \
	-lvtkCommon -lvtksys \
	-lvtkHybrid -lvtkexoIIc -lvtkNetCDF \
	-lvtkVolumeRendering -lvtkRendering -lvtkftgl \
	-lexpat $(VTK_EXPAT_LIB) $(VTK_FREETYPE_LIB) \
	-ltiff $(VTK_TIFF_LIB) -ljpeg $(VTK_JPEG_LIB) \
	$(VTK_PNG_LIB) $(VTK_Z_LIB) $(VTK_METAIO_LIB) \
	-lvtkIO	-lvtkDICOMParser $(VTK_MPEG2ENCODE_LIB) \
	-lvtkImagingTCL -lvtkVolumeRenderingTCL -lvtkRenderingTCL \
	-lvtkFilteringTCL -lvtkWidgetsTCL \
	-lvtkHybridTCL -lvtkGraphicsTCL -lvtkImagingTCL -lvtkIOTCL \
	-lvtkCommonTCL \
	$(KWWIDGETS_LIBS) -lKWWidgets

# vtkArrowPipelineTester test program
vtkArrowPipelineTester_SOURCES=vtkArrowPipelineTester.cxx

vtkArrowPipelineTester_LDADD=\
	$(top_builddir)/vtkutils/libvtkutils.a \
	$(LDFLAGS_TCL) $(LIBS_TCL_OPENGL) \
	$(VTK_LIBS) $(VTK_VERDICT_LIB) -lvtkWidgets -lvtkGraphics \
	-lvtkImaging -lvtkFiltering \
	-lvtkCommon -lvtksys \
	-lvtkHybrid -lvtkexoIIc -lvtkNetCDF \
	-lvtkVolumeRendering -lvtkRendering -lvtkftgl \
	-lexpat $(VTK_EXPAT_LIB) $(VTK_FREETYPE_LIB) \
	-ltiff $(VTK_TIFF_LIB) -ljpeg $(VTK_JPEG_LIB) $(VTK_SQLITE_LIB) \
	$(VTK_PNG_LIB) $(VTK_Z_LIB) $(VTK_METAIO_LIB) \
	-lvtkIO	-lvtkDICOMParser $(VTK_MPEG2ENCODE_LIB) \
	-lvtkImagingTCL -lvtkVolumeRenderingTCL -lvtkRenderingTCL \
	-lvtkFilteringTCL -lvtkWidgetsTCL \
	-lvtkHybridTCL -lvtkGraphicsTCL -lvtkImagingTCL -lvtkIOTCL \
	-lvtkCommonTCL 

clean-local:
	rm -f vtkKWBltGraphTcl.cxx vtkKWRGBATransferFunctionEditorTcl.cxx

EXTRA_DIST=IconLoader.h \
vtkArrowPipeline.h \
vtkFDTensorGlyph.h \
vtkFreesurferLookupTable.h \
vtkInflatePolyData.h \
vtkKWBltGraph.h \
vtkKWProgressDialog.h \
vtkKWRGBATransferFunctionEditor.h \
vtkKWRGBATransferFunctionEditorTester.h \
vtkODFGlyph.h \
vtkRGBATransferFunction.h \
vtkRGBAColorTransferFunction.h \
RgbatransferfunctioneditortesterLibInit.data \
vtkKWRGBATransferFunctionEditorTester-scalars.mgh \
IconLoaderTest.cxx IconLoaderTestIcons.txt testIconLoader

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra

clean-local:
	-rm *Tcl.cxx
