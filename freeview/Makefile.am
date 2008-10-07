AM_CXXFLAGS=\
	-fno-strict-aliasing \
	-Wno-deprecated \
	-I$(top_srcdir)/include \
	$(TCL_CFLAGS) \
	$(VTK_CXXFLAGS) \
	-I$(top_srcdir)/vtkfsio \
	-I$(top_srcdir)/vtkutils \
	$(WXWIDGETS_CXXFLAGS)

if ENABLE_WXWIDGETS_APPS
bin_PROGRAMS = freeview
freeview_SOURCES=\
	Annotation2D.cpp \
	Annotation2D.h \
	Broadcaster.cpp \
	Broadcaster.h \
	BrushProperty.cpp \
	BrushProperty.h \
	CommonDataStruct.h \
	ControlPanel.cpp \
	ControlPanel.h \
	Cursor2D.cpp \
	Cursor2D.h \
	Cursor3D.cpp \
	Cursor3D.h \
	CursorFactory.cpp \
	CursorFactory.h \
	DialogEditLookupTable.cpp \
	DialogEditLookupTable.h \
	DialogLoadDTI.cpp \
	DialogLoadDTI.h \
	DialogLoadVolume.cpp \
	DialogLoadVolume.h \
	DialogNewROI.cpp \
	DialogNewROI.h \
	DialogNewVolume.cpp \
	DialogNewVolume.h \
	DialogNewWayPoints.cpp \
	DialogNewWayPoints.h \
	DialogPreferences.cpp \
	DialogPreferences.h \
	DialogRotateVolume.cpp \
	DialogRotateVolume.h \
	FSLabel.cpp \
	FSLabel.h \
	FSSurface.cpp \
	FSSurface.h \
	FSVolume.cpp \
	FSVolume.h \
	FSWayPoints.cpp \
	FSWayPoints.h \
	Interactor.cpp \
	Interactor.h \
	Interactor2D.cpp \
	Interactor2D.h \
	Interactor2DNavigate.cpp \
	Interactor2DNavigate.h \
	Interactor2DRegionEdit.cpp \
	Interactor2DRegionEdit.h \
	Interactor2DROIEdit.cpp \
	Interactor2DROIEdit.h \
	Interactor2DVoxelEdit.cpp \
	Interactor2DVoxelEdit.h \
	Interactor2DWayPointsEdit.cpp \
	Interactor2DWayPointsEdit.h \
	Interactor3D.cpp \
	Interactor3D.h \
	Interactor3DNavigate.cpp \
	Interactor3DNavigate.h \
	LayerCollection.cpp \
	LayerCollection.h \
	LayerCollectionManager.cpp \
	LayerCollectionManager.h \
	Layer.cpp \
	Layer.h \
	LayerDTI.cpp \
	LayerDTI.h \
	LayerEditable.cpp \
	LayerEditable.h \
	LayerMRI.cpp \
	LayerMRI.h \
	LayerPropertiesDTI.cpp \
	LayerPropertiesDTI.h \
	LayerPropertiesMRI.cpp \
	LayerPropertiesMRI.h \
	LayerPropertiesROI.cpp \
	LayerPropertiesROI.h \
	LayerPropertiesSurface.cpp \
	LayerPropertiesSurface.h \
	LayerPropertiesWayPoints.cpp \
	LayerPropertiesWayPoints.h \
	LayerROI.cpp \
	LayerROI.h \
	LayerSurface.cpp \
	LayerSurface.h \
	LayerVolumeBase.cpp \
	LayerVolumeBase.h \
	LayerWayPoints.cpp \
	LayerWayPoints.h \
	Listener.cpp \
	Listener.h \
	LUTDataHolder.cpp \
	LUTDataHolder.h \
	MainApp.cpp \
	MainApp.h \
	MainWindow.cpp \
	MainWindow.h \
	MyUtils.cpp \
	MyUtils.h \
	MyCmdLineParser.cpp \
	MyCmdLineParser.h \
	PanelROI.cpp \
	PanelROI.h \
	PanelSceneSetting.cpp \
	PanelSceneSetting.h \
	PanelSurface.cpp \
	PanelSurface.h \
	PanelVolume.cpp \
	PanelVolume.h \
	PanelWayPoints.cpp \
	PanelWayPoints.h \
	PixelInfoPanel.cpp \
	PixelInfoPanel.h \
	Region.cpp \
	Region.h \
	RenderView2D.cpp \
	RenderView2D.h \
	RenderView3D.cpp \
	RenderView3D.h \
	RenderView.cpp \
	RenderView.h \
	StatusBar.cpp \
	StatusBar.h \
	string_fixed.h \
	ToolWindowEdit.cpp \
	ToolWindowEdit.h \
	WindowQuickReference.cpp \
	WindowQuickReference.h \
	WorkerThread.cpp \
	WorkerThread.h \
	wxColorIndicator.cpp \
	wxColorIndicator.h \
	wxColorIndicatorXmlHandler.cpp \
	wxColorIndicatorXmlHandler.h \
	wxVTKRenderWindowInteractor.cpp \
	wxVTKRenderWindowInteractor.h \
	res/resource.cpp \
	res/resource.h \
	res/QuickRef.h \
	res/CursorFill_png.h \
	res/CursorPencil_png.h \
	res/CursorPolyline_png.h \
	res/CursorPan_png.h \
	res/CursorZoom_png.h

freeview_LDFLAGS=$(OS_LDFLAGS) $(LDFLAGS_TCL)

freeview_LDADD=\
	$(addprefix $(top_builddir)/, $(LIBS_MGH)) $(LIBS_TCL_OPENGL) \
	$(top_builddir)/vtkfsio/libvtkfsio.a \
	$(top_builddir)/vtkutils/libvtkutils.a \
	$(VTK_LIBS) $(VTK_VERDICT_LIB) -lvtkGraphics \
	$(VTK_EXPAT_LIB) $(VTK_FREETYPE_LIB) \
	$(VTK_TIFF_LIB) $(VTK_JPEG_LIB) \
	$(VTK_METAIO_LIB) $(VTK_PNG_LIB) $(VTK_Z_LIB) \
	$(VTK_SQLITE_LIB) \
	-lvtkImaging -lvtkFiltering \
	-lvtkCommon -lvtksys -lvtkGenericFiltering \
	-lvtkexoIIc -lvtkNetCDF \
	-lvtkVolumeRendering -lvtkRendering -lvtkftgl \
	-lvtkWidgets -lvtkHybrid \
	-lvtkIO	-lvtkDICOMParser $(VTK_MPEG2ENCODE_LIB) \
	$(WXWIDGETS_LIBS)

# put a wrapper around the bin, used to setup tcltktixblt,vtk,kww enviro vars
install-exec-hook:freeview
	cp $(top_builddir)/freeview/freeview $(DESTDIR)$(bindir)/freeview.bin
	echo "#!/bin/tcsh -ef" > $(DESTDIR)$(bindir)/freeview
	echo "source \$$FREESURFER_HOME/bin/tcl_setup" \
		>> $(DESTDIR)$(bindir)/freeview
	echo "source \$$FREESURFER_HOME/bin/vtk_setup" \
		>> $(DESTDIR)$(bindir)/freeview
	echo "freeview.bin \$$argv" >> $(DESTDIR)$(bindir)/freeview
	chmod a+x $(DESTDIR)$(bindir)/freeview
	if [[ -f /usr/pubsw/bin/chrpath ]]; then \
	/usr/pubsw/bin/chrpath -d $(DESTDIR)$(bindir)/freeview.bin; fi

uninstall-hook:
	rm -f $(DESTDIR)$(bindir)/freeview.bin

else
bin_PROGRAMS =

all:
	@echo "Building wxWidgets apps not enabled."
	@echo "Run configure with either --enable-wxw-apps or --with-wxwidgets-dir=DIR"

endif

EXTRA_DIST=\
	res/DialogEditLookupTable.xrc \
	res/DialogLoadDTI.xrc \
	res/DialogLoadVolume.xrc \
	res/DialogNewROI.xrc \
	res/DialogNewVolume.xrc \
	res/DialogPreferences.xrc \
	res/MainWindow.xrc \
	res/PanelROI.xrc \
	res/PanelVolume.xrc \
	res/WindowQuickReference.xrc \
  res/QuickRef.html \
	res/images/edit_redo_16.png \
	res/images/edit_redo.png \
	res/images/edit_undo_16.png \
	res/images/edit_undo.png \
	res/images/file_open_16.png \
	res/images/file_open.png \
	res/images/layout_1n3.png \
	res/images/layout_1x1.png \
	res/images/layout_2x2.png \
	res/images/move_down.png \
	res/images/move_up.png \
	res/images/navigate.png \
	res/images/roi_delete.png \
	res/images/roi_draw_fill.png \
	res/images/roi_draw_freehand.png \
	res/images/roi_draw_polyline.png \
	res/images/roi_edit.png \
	res/images/roi_load.png \
	res/images/roi_new.png \
	res/images/roi_save.png \
	res/images/view_3d.png \
	res/images/view_axial.png \
	res/images/view_coronal.png \
	res/images/view_sagittal.png \
	res/images/volume_delete.png \
	res/images/volume_load_16.png \
	res/images/volume_load.png \
	res/images/volume_new.png \
	res/images/volume_save.png \
	res/images/voxel_draw_fill.png \
	res/images/voxel_draw_freehand.png \
	res/images/voxel_draw_polyline.png \
	res/images/voxel_edit.png

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
