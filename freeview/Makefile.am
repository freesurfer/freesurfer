# The logic is divided into 2 parts based on whether we are building for OSX or not.
# The VTK for OSX doesn't compile when Carbon is enabled and TCL is enabled. However,
# a TCL enabled VTK is essential for FreeSurfer suite. Hence Freeview uses a special
# static build of VTK which has Carbon=ON, TCL=OFF, Static=ON. So depending on whether
# the platform is OSX or not, we use different versions of the VTK. 

if ENABLE_WXWIDGETS_APPS
bin_PROGRAMS = freeview
freeview_SOURCES=\
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

freeview_LDFLAGS=$(OS_LDFLAGS) $(LDFLAGS_TCL)

# put a wrapper around the bin, used to setup tcltktixblt,vtk,kww enviro vars
# and if the OS is MACOSX have a different install-hook to create a bundle
freeview_LDADD=\
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

install-exec-hook:freeview
	mkdir -p $(DESTDIR)$(bindir)/Freeview.app
	mkdir -p $(DESTDIR)$(bindir)/Freeview.app/Contents
	mkdir -p $(DESTDIR)$(bindir)/Freeview.app/Contents/MacOS
	mkdir -p $(DESTDIR)$(bindir)/Freeview.app/Contents/Resources
	mkdir -p $(DESTDIR)$(bindir)/Freeview.app/Contents/Resources/English.lproj
#	mkdir -p $(DESTDIR)$(bindir)/Freeview.app/Contents/Frameworks
	cp $(top_builddir)/freeview/Info.plist $(DESTDIR)$(bindir)/Freeview.app/Contents
	echo -n 'APPL????' > $(DESTDIR)$(bindir)/Freeview.app/Contents/PkgInfo
	cp $(top_builddir)/freeview/freeview $(DESTDIR)$(bindir)/Freeview.app/Contents/MacOS/freeview.bin
	echo "#!/bin/tcsh -ef" > $(DESTDIR)$(bindir)/freeview
#	echo "source \$$FREESURFER_HOME/bin/tcl_setup" \
		>> $(DESTDIR)$(bindir)/freeview
#	echo "source \$$FREESURFER_HOME/bin/vtk_setup" \
		>> $(DESTDIR)$(bindir)/freeview
	echo "exec $(DESTDIR)$(bindir)/Freeview.app/Contents/MacOS/freeview.bin \$$*" >> $(DESTDIR)$(bindir)/freeview
	chmod a+x $(DESTDIR)$(bindir)/freeview
#	./mac_package_dylibs.sh $(WXWIDGETS_DIR)lib $(DESTDIR)$(bindir)/Freeview.app/Contents/Frameworks $(DESTDIR)$(bindir)/Freeview.app/Contents/MacOS/freeview.bin
	if [[ -f /usr/pubsw/bin/chrpath ]]; then \
	/usr/pubsw/bin/chrpath -d $(DESTDIR)$(bindir)/Freeview.app/Contents/MacOS/freeview.bin; fi

uninstall-hook:
	rm -f $(DESTDIR)$(bindir)/freeview
	rm -rf $(DESTDIR)$(bindir)/Freeview.app
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

freeview_LDFLAGS=$(OS_LDFLAGS) $(LDFLAGS_TCL)

freeview_LDADD=\
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
endif

else
bin_PROGRAMS =

all:
	@echo "Building wxWidgets apps not enabled."
	@echo "Run configure with either --enable-wxw-apps or --with-wxwidgets-dir=DIR"

endif

# NJS commented-out till a virtual-screen way of testing is implemented
#TESTS=test_freeview

EXTRA_DIST=test_freeview

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
