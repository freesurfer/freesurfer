#
# Makefile.am for libvtkfsio.a
#

AM_CPPFLAGS= -Wno-deprecated \
	$(X_CFLAGS) \
	$(VXL_CFLAGS) \
	-I$(top_srcdir)/include \
	$(VTK_CXXFLAGS) \
	$(KWWIDGETS_CXXFLAGS)

AM_LDFLAGS=

noinst_LIBRARIES= libvtkfsio.a

libvtkfsio_a_SOURCES=\
	vtkFSSurfaceLabelSource.cxx \
	vtkFSSurfaceWriter.cxx \
	vtkFSSurfaceScalarsReader.cxx \
	vtkFSSurfaceSource.cxx \
	vtkFSVolumeSource.cxx

am_libvtkfsio_a_OBJECTS=$(libvtkfsio_a_SOURCES:%.cxx=%.o)

EXTRA_DIST=$(wildcard *.h)

if ENABLE_KWWIDGETS_APPS
check_PROGRAMS = vtkFSSurfaceWriterTester
else
check_PROGRAMS =
endif

vtkFSSurfaceWriterTester_SOURCES=vtkFSSurfaceWriterTester.cxx

vtkFSSurfaceWriterTester_LDADD=\
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
	$(VTK_PNG_LIB) $(VTK_Z_LIB) $(VTK_METAIO_LIB) \
	-lvtkIO	-lvtkDICOMParser $(VTK_MPEG2ENCODE_LIB) \
	-lvtkImagingTCL -lvtkVolumeRenderingTCL -lvtkRenderingTCL \
	-lvtkFilteringTCL -lvtkWidgetsTCL \
	-lvtkHybridTCL -lvtkGraphicsTCL -lvtkImagingTCL -lvtkIOTCL \
	-lvtkCommonTCL

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
