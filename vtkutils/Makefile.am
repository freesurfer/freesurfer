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
	vtkKWProgressDialog.cxx vtkKWProgressDialog.h \
	vtkKWRGBATransferFunctionEditor.cxx vtkKWRGBATransferFunctionEditor.h \
	vtkKWRGBATransferFunctionEditorTcl.cxx \
	vtkFreesurferLookupTable.cxx vtkFreesurferLookupTable.h \
	vtkRGBATransferFunction.cxx vtkRGBATransferFunction.h \
	vtkFDTensorGlyph.cxx vtkFDTensorGlyph.h \
	vtkInflatePolyData.cxx vtkInflatePolyData.h \
	IconLoader.cxx IconLoader.h

am_libvtkutils_a_OBJECTS=$(libvtkutils_a_SOURCES:%.cxx=%.o)

# We use vtkWrapTcl to create the Tcl wrapped code for this object.
vtkKWRGBATransferFunctionEditorTcl.cxx: vtkKWRGBATransferFunctionEditor.h
	${VTK_DIR}/bin/vtkWrapTcl vtkKWRGBATransferFunctionEditor.h \
	${VTK_LIB_VTK}/hints 1 vtkKWRGBATransferFunctionEditorTcl.cxx


if ENABLE_KWWIDGETS_APPS
check_PROGRAMS=vtkKWRGBATransferFunctionEditorTester
else
check_PROGRAMS=
endif

#TESTS=$(check_PROGRAMS)

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
	$(top_srcdir)/fsgdf/libfsgdf.a \
	$(top_srcdir)/vtkfsio/libvtkfsio.a \
	$(top_srcdir)/vtkutils/libvtkutils.a \
	$(addprefix $(top_srcdir)/, $(LIBS_MGH)) $(LDFLAGS_TCL) $(LIBS_TCL_OPENGL) \
	$(VTK_LIBS) $(VTK_VERDICT_LIB) -lvtkWidgets -lvtkGraphics \
	-lvtkImaging -lvtkFiltering \
	-lvtkCommon -lvtksys \
	-lvtkHybrid -lvtkexoIIc -lvtkNetCDF \
	-lvtkVolumeRendering -lvtkRendering -lvtkftgl \
	$(VTK_EXPAT_LIB) $(VTK_FREETYPE_LIB) \
	$(VTK_TIFF_LIB) $(VTK_JPEG_LIB) \
	$(VTK_PNG_LIB) $(VTK_Z_LIB) \
	-lvtkIO	-lvtkDICOMParser $(VTK_MPEG2ENCODE_LIB) \
	-lvtkImagingTCL -lvtkVolumeRenderingTCL -lvtkRenderingTCL \
	-lvtkFilteringTCL -lvtkWidgetsTCL \
	-lvtkHybridTCL -lvtkGraphicsTCL -lvtkImagingTCL -lvtkIOTCL \
	-lvtkCommonTCL \
	$(KWWIDGETS_LIBS) -lKWWidgets

EXTRA_DIST=$(wildcard *.h) $(wildcard *.data)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
