#
# Makefile.am for libvtkutils.a
#

AM_CPPFLAGS= -Wno-deprecated \
	$(X_CFLAGS) \
	$(TCL_CFLAGS) \
	$(VXL_CFLAGS) \
	-I$(top_srcdir)/include \
	$(VTK_CXXFLAGS) \
	$(KWWIDGETS_CXXFLAGS)

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

check_PROGRAMS=foo vtkKWRGBATransferFunctionEditorTester vtkArrowPipelineTester
#TESTS is commented-out because the tests are interactive (can't run automatically)
#TESTS=vtkKWRGBATransferFunctionEditorTester vtkArrowPipelineTester
else
check_PROGRAMS=
TESTS=
endif

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
	$(addprefix $(top_builddir)/, $(LIBS_MGH)) $(LDFLAGS_TCL) $(LIBS_TCL_OPENGL) \
	$(VTK_LIBS) $(VTK_VERDICT_LIB) -lvtkWidgets -lvtkGraphics \
	-lvtkImaging -lvtkFiltering \
	-lvtkCommon -lvtksys \
	-lvtkHybrid -lvtkexoIIc -lvtkNetCDF \
	-lvtkVolumeRendering -lvtkRendering -lvtkftgl \
	$(VTK_EXPAT_LIB) $(VTK_FREETYPE_LIB) \
	$(VTK_TIFF_LIB) $(VTK_JPEG_LIB) $(VTL_SQLITE_LIB) \
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
	$(VTK_LIBS) $(VTK_VERDICT_LIB) -lvtkWidgets -lvtkGraphics \
	-lvtkImaging -lvtkFiltering \
	-lvtkCommon -lvtksys \
	-lvtkHybrid -lvtkexoIIc -lvtkNetCDF \
	-lvtkVolumeRendering -lvtkRendering -lvtkftgl \
	$(VTK_EXPAT_LIB) $(VTK_FREETYPE_LIB) \
	$(VTK_TIFF_LIB) $(VTK_JPEG_LIB) $(VTL_SQLITE_LIB) \
	$(VTK_PNG_LIB) $(VTK_Z_LIB) $(VTK_METAIO_LIB) \
	-lvtkIO	-lvtkDICOMParser $(VTK_MPEG2ENCODE_LIB) \
	-lvtkImagingTCL -lvtkVolumeRenderingTCL -lvtkRenderingTCL \
	-lvtkFilteringTCL -lvtkWidgetsTCL \
	-lvtkHybridTCL -lvtkGraphicsTCL -lvtkImagingTCL -lvtkIOTCL \
	-lvtkCommonTCL 


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
RgbatransferfunctioneditortesterLibInit.data \
vtkKWRGBATransferFunctionEditorTester-scalars.mgh

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
