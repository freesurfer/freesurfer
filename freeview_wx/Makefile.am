# The logic is divided into 2 parts based on whether we are building for 
# OSX or not. The VTK for OSX doesn't compile when Carbon is enabled and 
# TCL is enabled. However, a TCL enabled VTK is essential for FreeSurfer
# suite. Hence Freeview uses a special static build of VTK which has 
# Carbon=ON, TCL=OFF, Static=ON. So depending on whether the platform is 
# OSX or not, we use different versions of the VTK. 

if ENABLE_WXWIDGETS_APPS
bin_PROGRAMS = freeview_wx
freeview_wx_SOURCES=\
	Annotation2D.cpp \
	Annotation2D.h \
	BuildContourThread.cpp \
	BuildContourThread.h \
	Broadcaster.cpp \
	Broadcaster.h \
	BrushProperty.cpp \
	BrushProperty.h \
	CommonDataStruct.h \
	ConnectivityData.cpp \
	ConnectivityData.h \
	ControlPanel.cpp \
	ControlPanel.h \
	Contour2D.cpp \
	Contour2D.h \
	Cursor2D.cpp \
	Cursor2D.h \
	Cursor3D.cpp \
	Cursor3D.h \
	CursorFactory.cpp \
	CursorFactory.h \
	DialogCropVolume.cpp \
	DialogCropVolume.h \
	DialogEditLookupTable.cpp \
	DialogEditLookupTable.h \
	DialogGradientVolume.cpp \
	DialogGradientVolume.h \
	DialogLoadDTI.cpp \
	DialogLoadDTI.h \
	DialogLoadROI.cpp \
	DialogLoadROI.h \
	DialogLoadPointSet.cpp \
	DialogLoadPointSet.h \
	DialogLoadPVolumes.cpp \
	DialogLoadPVolumes.h \
	DialogLoadVolume.cpp \
	DialogLoadVolume.h \
	DialogNewROI.cpp \
	DialogNewROI.h \
	DialogNewVolume.cpp \
	DialogNewVolume.h \
	DialogNewWayPoints.cpp \
	DialogNewWayPoints.h \
	DialogOptimalVolume.cpp \
	DialogOptimalVolume.h \
	DialogPreferences.cpp \
	DialogPreferences.h \
	DialogRepositionSurface.cpp \
	DialogRepositionSurface.h \
	DialogSavePoint.cpp \
	DialogSavePoint.h \
	DialogSavePointSetAs.cpp \
	DialogSavePointSetAs.h \
	DialogSaveScreenshot.cpp \
	DialogSaveScreenshot.h \
	DialogSaveVolumeAs.cpp \
	DialogSaveVolumeAs.h \
	DialogTransformVolume.cpp \
	DialogTransformVolume.h \
	DialogVolumeFilter.cpp \
	DialogVolumeFilter.h \
	DialogWriteMovieFrames.cpp \
	DialogWriteMovieFrames.h \
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
	Interactor2DCropVolume.cpp \
	Interactor2DCropVolume.h \
	Interactor2DNavigate.cpp \
	Interactor2DNavigate.h \
	Interactor2DMeasure.cpp \
	Interactor2DMeasure.h \
	Interactor2DRegionEdit.cpp \
	Interactor2DRegionEdit.h \
	Interactor2DROIEdit.cpp \
	Interactor2DROIEdit.h \
	Interactor2DVolumeEdit.cpp \
	Interactor2DVolumeEdit.h \
	Interactor2DVoxelEdit.cpp \
	Interactor2DVoxelEdit.h \
	Interactor2DWayPointsEdit.cpp \
	Interactor2DWayPointsEdit.h \
	Interactor3D.cpp \
	Interactor3D.h \
	Interactor3DCropVolume.cpp \
	Interactor3DCropVolume.h \
	Interactor3DMeasure.cpp \
	Interactor3DMeasure.h \
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
	LayerOptimal.cpp \
	LayerOptimal.h \
	LayerPLabel.cpp \
	LayerPLabel.h \
	LayerProperties.cpp \
	LayerProperties.h \
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
	LivewireTool.cpp \
	LivewireTool.h \
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
	Region2D.cpp \
	Region2D.h \
	Region2DLine.cpp \
	Region2DLine.h \
	Region2DPolyline.cpp \
	Region2DPolyline.h \
	Region2DRectangle.cpp \
	Region2DRectangle.h \
	RenderView2D.cpp \
	RenderView2D.h \
	RenderView3D.cpp \
	RenderView3D.h \
	RenderView.cpp \
	RenderView.h \
	StatusBar.cpp \
	StatusBar.h \
	string_fixed.h \
	SurfaceAnnotation.h \
	SurfaceAnnotation.cpp \
	SurfaceLabel.h \
	SurfaceLabel.cpp \
	SurfaceOverlay.cpp \
	SurfaceOverlay.h \
	SurfaceOverlayProperties.cpp \
	SurfaceOverlayProperties.h \
	SurfaceRegion.cpp \
	SurfaceRegion.h \
	SurfaceRegionGroups.cpp \
	SurfaceRegionGroups.h \
	ToolWindowEdit.cpp \
	ToolWindowEdit.h \
	ToolWindowMeasure.cpp \
	ToolWindowMeasure.h \
	VolumeCropper.cpp \
	VolumeCropper.h \
	VolumeFilter.cpp \
	VolumeFilter.h \
	VolumeFilterConvolve.cpp \
	VolumeFilterConvolve.h \
	VolumeFilterGradient.cpp \
	VolumeFilterGradient.h \
	VolumeFilterMean.cpp \
	VolumeFilterMean.h \
	VolumeFilterMedian.cpp \
	VolumeFilterMedian.h \
	vtkSimpleLabelEdgeFilter.cpp \
	vtkSimpleLabelEdgeFilter.h \
	WindowConnectivityConfiguration.cpp \
	WindowConnectivityConfiguration.h \
	WindowHistogram.cpp \
	WindowHistogram.h \
	WindowQuickReference.cpp \
	WindowQuickReference.h \
	WindowOverlayConfiguration.cpp \
	WindowOverlayConfiguration.h \
	WorkerThread.cpp \
	WorkerThread.h \
	wxColorIndicator.cpp \
	wxColorIndicator.h \
	wxColorIndicatorXmlHandler.cpp \
	wxColorIndicatorXmlHandler.h \
	wxHistogramWidget.cpp \
	wxHistogramWidget.h \
	wxHistogramWidgetXmlHandler.cpp \
	wxHistogramWidgetXmlHandler.h \
	wxVTKRenderWindowInteractor.cpp \
	wxVTKRenderWindowInteractor.h \
	res/resource.cpp \
	res/resource.h \
	res/QuickRef.h

RESOURCES=\
	res/DialogLoadDTI.xrc \
	res/DialogVolumeFilter.xrc \
	res/CursorFill_png.h \
	res/CursorAdd_png.h \
	res/DialogNewVolume.xrc \
	res/DialogPreferences.xrc \
	res/DialogLoadROI.xrc \
	res/DialogRepositionSurface.xrc \
	res/DialogSavePoint.xrc \
	res/PanelSurface.xrc \
	res/CursorZoom_png.h \
	res/QuickRef.html \
	res/QuickRef.h \
	res/DialogOptimalVolume.xrc \
	res/DialogSavePointSetAs.xrc \
	res/DialogEditLookupTable.xrc \
	res/ToolWindowMeasure.xrc \
	res/ToolWindowEdit.xrc \
	res/PanelVolume.xrc \
	res/FreeSurferColorLUT.h \
	res/DialogTransformVolume.xrc \
	res/CursorMeasure_png.h \
	res/resource.cpp \
	res/DialogLoadConnectivityData.xrc \
	res/CursorPolyline_png.h \
	res/DialogSaveVolumeAs.xrc \
	res/CursorMeasurePolyline_png.h \
	res/DialogGradientVolume.xrc \
	res/DialogSaveScreenshot.xrc \
	res/CursorMeasureLine_png.h \
	res/CursorPan_png.h \
	res/DialogNewROI.xrc \
	res/res2c.sh \
	res/DialogLoadPVolumes.xrc \
	res/WindowQuickReference.xrc \
	res/CursorColorPicker_png.h \
	res/CursorGrab_png.h \
	res/CursorRemove_png.h \
	res/PanelSceneSetting.xrc \
	res/MainWindow.xrc \
	res/DialogWriteMovieFrames.xrc \
	res/resource.h \
	res/PanelROI.xrc \
	res/DialogLoadPointSet.xrc \
	res/CursorPencil_png.h \
	res/DialogNewWayPoints.xrc \
	res/images/surface_inflated.png \
	res/images/tool_save_point.png \
	res/images/view_histogram.png \
	res/images/surface_main.png \
	res/images/edit_undo.png \
	res/images/surface_pial.png \
	res/images/voxel_draw_freehand.png \
	res/images/view_control_panel.png \
	res/images/surface_load_16.png \
	res/images/cursor_zoom.png \
	res/images/waypoints_load.png \
	res/images/surface_load.png \
	res/images/cursor_fill.png \
	res/images/view_coordinate.png \
	res/images/view_sagittal.png \
	res/images/cursor_measure_polyline.png \
	res/images/voxel_draw_livewire.png \
	res/images/move_up.png \
	res/images/measure_line.png \
	res/images/roi_save.png \
	res/images/edit_paste_16.png \
	res/images/waypoints_save.png \
	res/images/edit_redo_16.png \
	res/images/surface_curvature.png \
	res/images/view_scalar_bar.png \
	res/images/layout_1n3_h.png \
	res/images/measure_label.png \
	res/images/tool_goto_point.png \
	res/images/navigate.png \
	res/images/volume_paste_setting.png \
	res/images/cursor_polyline.png \
	res/images/volume_lock.png \
	res/images/edit_copy_16.png \
	res/images/view_3d.png \
	res/images/surface_original.png \
	res/images/file_open.png \
	res/images/volume_paste_setting_all.png \
	res/images/view_coronal.png \
	res/images/view_coronal_old.png \
	res/images/edit_redo.png \
	res/images/voxel_draw_fill.png \
	res/images/roi_load.png \
	res/images/measure_polyline.png \
	res/images/volume_save.png \
	res/images/move_down.png \
	res/images/surface_delete.png \
	res/images/cursor_measure_rectangle.png \
	res/images/view_3d_old.png \
	res/images/view_axial_old.png \
	res/images/cursor_remove.png \
	res/images/view_axial.png \
	res/images/volume_delete.png \
	res/images/layout_2x2.png \
	res/images/measure_rectangle.png \
	res/images/layout_1n3.png \
	res/images/edit_paste.png \
	res/images/measure.png \
	res/images/waypoints_new.png \
	res/images/view_sagittal_old.png \
	res/images/roi_draw_freehand.png \
	res/images/voxel_edit.png \
	res/images/cursor_colorpicker.png \
	res/images/measure_spline.png \
	res/images/roi_delete.png \
	res/images/layout_1x1.png \
	res/images/volume_new.png \
	res/images/file_screenshot_16.png \
	res/images/cursor_pan.png \
	res/images/cursor_contour.png \
	res/images/voxel_draw_colorpicker.png \
	res/images/surface_white.png \
	res/images/cursor_grab.png \
	res/images/voxel_draw_contour.png \
	res/images/cursor_measure_line.png \
	res/images/roi_draw_polyline.png \
	res/images/voxel_draw_polyline.png \
	res/images/file_open_16.png \
	res/images/volume_load.png \
	res/images/edit_copy.png \
	res/images/tool_rotate_volume.png \
	res/images/file_screenshot.png \
	res/images/measure_surface_region.png \
	res/images/roi_new.png \
	res/images/cursor_pencil.png \
	res/images/volume_load_16.png \
	res/images/roi_draw_fill.png \
	res/images/waypoints_edit.png \
	res/images/volume_copy_setting.png \
	res/images/edit_undo_16.png \
	res/images/roi_edit.png \
	res/images/cursor_add.png \
	res/images/waypoints_delete.png \
	res/PanelWayPoints.xrc \
	res/DialogLoadVolume.xrc \
	res/CursorMeasureRectangle_png.h \
	res/CursorContour_png.h \
	res/WindowOverlayConfiguration.xrc \
	res/WindowHistogram.xrc \
	res/WindowConnectivityConfiguration.xrc \
	res/Info.plist \
	res/DialogCropVolume.xrc

if HAVE_MAC_OSX
AM_CXXFLAGS=\
	-fno-strict-aliasing \
	-Wno-deprecated \
	-I$(top_srcdir)/include \
	$(TCL_CFLAGS) \
	$(GL_CFLAGS) \
	$(VTK_CARBON_CXXFLAGS) \
	-I$(top_srcdir)/vtkfsio \
	-I$(top_srcdir)/vtkutils \
	$(WXWIDGETS_CXXFLAGS) \
	$(ITK_CFLAGS) -DHAVE_ITK_LIBS \
	$(VXL_CFLAGS)

freeview_wx_LDFLAGS=$(OS_LDFLAGS) $(LDFLAGS_TCL)

# put a wrapper around the bin, used to setup tcltktixblt,vtk,kww enviro vars
# and if the OS is MACOSX have a different install-hook to create a bundle
freeview_wx_LDADD=\
	$(addprefix $(top_builddir)/, $(LIBS_MGH)) \
	$(top_builddir)/vtkfsio/libvtkfsio.a \
	$(top_builddir)/vtkutils/libvtkutils.a \
	$(VTK_CARBON_LIBS) $(VTK_VERDICT_LIB) -lvtkGraphics \
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
	$(WXWIDGETS_LIBS) $(WXGTK2GL) $(WXWIDGETS_OPTIONAL_LIBS)

install-exec-hook:freeview_wx
	mkdir -p $(DESTDIR)$(bindir)/Freeview.wx.app
	mkdir -p $(DESTDIR)$(bindir)/Freeview.wx.app/Contents
	mkdir -p $(DESTDIR)$(bindir)/Freeview.wx.app/Contents/MacOS
	mkdir -p $(DESTDIR)$(bindir)/Freeview.wx.app/Contents/Resources
	mkdir -p $(DESTDIR)$(bindir)/Freeview.wx.app/Contents/Resources/English.lproj
#	mkdir -p $(DESTDIR)$(bindir)/Freeview.wx.app/Contents/Frameworks
	cp $(top_builddir)/freeview_wx/Info.plist $(DESTDIR)$(bindir)/Freeview.wx.app/Contents
	echo -n 'APPL????' > $(DESTDIR)$(bindir)/Freeview.wx.app/Contents/PkgInfo
	cp $(top_builddir)/freeview_wx/freeview_wx $(DESTDIR)$(bindir)/Freeview.wx.app/Contents/MacOS/freeview_wx.bin
	echo "#!/bin/tcsh -ef" > $(DESTDIR)$(bindir)/freeview_wx
#	echo "source \$$FREESURFER_HOME/bin/tcl_setup" \
		>> $(DESTDIR)$(bindir)/freeview_wx
#	echo "source \$$FREESURFER_HOME/bin/vtk_setup" \
		>> $(DESTDIR)$(bindir)/freeview_wx
	echo "exec $(DESTDIR)$(bindir)/Freeview.wx.app/Contents/MacOS/freeview_wx.bin \$$*" >> $(DESTDIR)$(bindir)/freeview_wx
	chmod a+x $(DESTDIR)$(bindir)/freeview_wx
#	./mac_package_dylibs.sh $(WXWIDGETS_DIR)lib $(DESTDIR)$(bindir)/Freeview.wx.app/Contents/Frameworks $(DESTDIR)$(bindir)/Freeview.wx.app/Contents/MacOS/freeview_wx.bin
	if [[ -f /usr/pubsw/bin/chrpath ]]; then \
	/usr/pubsw/bin/chrpath -d $(DESTDIR)$(bindir)/Freeview.wx.app/Contents/MacOS/freeview_wx.bin; fi

uninstall-hook:
	rm -f $(DESTDIR)$(bindir)/freeview_wx
	rm -rf $(DESTDIR)$(bindir)/Freeview.wx.app
else
AM_CXXFLAGS=\
	-fno-strict-aliasing \
	-Wno-deprecated \
	-I$(top_srcdir)/include \
	$(TCL_CFLAGS) \
	$(GL_CFLAGS) \
	$(VTK_CXXFLAGS) \
	-I$(top_srcdir)/vtkfsio \
	-I$(top_srcdir)/vtkutils \
	$(WXWIDGETS_CXXFLAGS) \
	$(ITK_CFLAGS) -DHAVE_ITK_LIBS \
	$(VXL_CFLAGS)

freeview_wx_LDFLAGS=$(OS_LDFLAGS) $(LDFLAGS_TCL)

freeview_wx_LDADD=\
	$(addprefix $(top_builddir)/, $(LIBS_MGH)) \
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
	$(WXWIDGETS_LIBS) $(WXGTK2GL) $(WXWIDGETS_OPTIONAL_LIBS)

# put a wrapper around the bin, used to setup tcltktixblt,vtk,kww enviro vars
# and if the OS is MACOSX have a different install-hook to create a bundle
install-exec-hook:freeview_wx
	cp $(top_builddir)/freeview_wx/freeview_wx $(DESTDIR)$(bindir)/freeview_wx.bin
	echo "#!/bin/tcsh -ef" > $(DESTDIR)$(bindir)/freeview_wx
	echo "source \$$FREESURFER_HOME/bin/tcl_setup" \
		>> $(DESTDIR)$(bindir)/freeview_wx
	echo "source \$$FREESURFER_HOME/bin/vtk_setup" \
		>> $(DESTDIR)$(bindir)/freeview_wx
	echo "freeview_wx.bin \$$argv" >> $(DESTDIR)$(bindir)/freeview_wx
	chmod a+x $(DESTDIR)$(bindir)/freeview_wx
	if [[ -f /usr/pubsw/bin/chrpath ]]; then \
	/usr/pubsw/bin/chrpath -d $(DESTDIR)$(bindir)/freeview_wx.bin; fi

uninstall-hook:
	rm -f $(DESTDIR)$(bindir)/freeview_wx.bin
endif

else
bin_PROGRAMS =

all:
	@echo "Building wxWidgets apps not enabled."
	@echo "Run configure with either --enable-wxw-apps or --with-wxwidgets-dir=DIR"

endif

# NJS commented-out till a virtual-screen way of testing is implemented
#TESTS=test_freeview_wx

EXTRA_DIST=test_freeview $(RESOURCES) StockColorMap.h Info.plist

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
