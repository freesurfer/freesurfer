
AM_CFLAGS=-I$(top_srcdir)/include
AM_CPPFLAGS=-I$(top_srcdir)/include $(TCL_CFLAGS) \
	-I/usr/include/libxml2
AM_CXXFLAGS=-Wno-deprecated -I$(top_srcdir)/include \
	-I$(top_srcdir)/vtkfsio \
	-I$(top_srcdir)/vtkutils \
	$(VTK_CXXFLAGS) $(KWWIDGETS_CXXFLAGS) $(TCL_CFLAGS)

bin_PROGRAMS = mris_mesh_subdivide
mris_mesh_subdivide_SOURCES=mris_mesh_subdivide.c
mris_mesh_subdivide_LDADD= \
	$(top_builddir)/vtkfsio/libvtkfsio.a \
	$(top_builddir)/vtkutils/libvtkutils.a \
	$(addprefix $(top_builddir)/, $(LIBS_MGH)) $(LIBS_TCL_OPENGL) \
	$(VTK_LIBS) $(VTK_VERDICT_LIB) -lvtkGraphics \
	$(VTK_EXPAT_LIB) $(VTK_FREETYPE_LIB) \
	$(VTK_TIFF_LIB) $(VTK_JPEG_LIB) \
	$(VTK_METAIO_LIB) $(VTK_PNG_LIB) $(VTK_Z_LIB) \
	$(VTK_SQLITE_LIB) \
	-lvtkImaging -lvtkFiltering \
	-lvtkCommon -lvtksys \
	-lvtkexoIIc -lvtkNetCDF \
	-lvtkVolumeRendering -lvtkRendering -lvtkftgl \
	-lvtkWidgets -lvtkHybrid \
	-lvtkIO	-lvtkDICOMParser $(VTK_MPEG2ENCODE_LIB) \
	-lvtkImagingTCL -lvtkVolumeRenderingTCL -lvtkRenderingTCL \
	-lvtkFilteringTCL -lvtkWidgetsTCL \
	-lvtkHybridTCL -lvtkGraphicsTCL -lvtkImagingTCL -lvtkIOTCL \
	-lvtkCommonTCL \
	$(KWWIDGETS_LIBS) -lKWWidgets
mris_mesh_subdivide_LDFLAGS=$(OS_LDFLAGS) $(LDFLAGS_TCL)



# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
