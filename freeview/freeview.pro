QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = freeview
TEMPLATE = app

SOURCES += \
    Annotation2D.cpp \
    CommandEdit.cpp \
    Contour2D.cpp \
    Cursor2D.cpp \
    Cursor3D.cpp \
    CursorFactory.cpp \
    BrushProperty.cpp  \
    DialogAbout.cpp \
    DialogControlPointComment.cpp \
    DialogCropVolume.cpp \
    DialogGradientFilter.cpp \
    DialogLoadDTI.cpp \
    DialogLoadODF.cpp \
    DialogLoadPointSet.cpp \
    DialogLoadVolume.cpp \
    DialogMovePoint.cpp \
    DialogNewAnnotation.cpp \
    DialogNewROI.cpp \
    DialogNewPointSet.cpp \
    DialogNewVolume.cpp \
    DialogPreferences.cpp \
    DialogSaveAllVolumes.cpp \
    DialogSavePointSet.cpp \
    DialogSaveScreenshot.cpp \
    DialogTransformSurface.cpp \
    DialogTransformVolume.cpp \
    DialogVolumeFilter.cpp \
    DialogWriteMovieFrames.cpp \
    FloatingStatusBar.cpp \
    FSLabel.cpp\
    FSPointSet.cpp \
    FSSurface.cpp \
    FSVolume.cpp \
    GenericRenderView.cpp \
    InfoTreeWidget.cpp \
    Interactor.cpp \
    Interactor2D.cpp \
    Interactor2DMeasure.cpp \
    Interactor2DNavigate.cpp \
    Interactor2DPointSetEdit.cpp \
    Interactor2DROIEdit.cpp \
    Interactor2DVolumeCrop.cpp \
    Interactor2DVolumeEdit.cpp \
    Interactor2DVoxelEdit.cpp \
    Interactor3D.cpp \
    Interactor3DMeasure.cpp \
    Interactor3DNavigate.cpp \
    Interactor3DVolumeCrop.cpp \
    Interactor3DPointSetEdit.cpp \
    Layer.cpp \
    LayerCollection.cpp \
    LayerDTI.cpp \
    LayerEditable.cpp \
    LayerMRI.cpp \
    LayerODF.cpp \
    LayerPLabel.cpp \
    LayerPointSet.cpp \
    LayerProperty.cpp \
    LayerPropertyDTI.cpp \
    LayerPropertyMRI.cpp \
    LayerPropertyODF.cpp \
    LayerPropertyPointSet.cpp \
    LayerPropertyROI.cpp \
    LayerPropertySurface.cpp \
    LayerROI.cpp \
    LayerSurface.cpp \
    LayerTreeWidget.cpp \
    LayerVolumeBase.cpp \
    LivewireTool.cpp \
    LUTDataHolder.cpp \
    PanelODF.cpp \
    Region3D.cpp \
    ScribblePromptWorker.cpp \
    ToolWindowLesionPopup.cpp \
    TorchScriptModule.cpp \
    VolumeFilterOptimal.cpp \
    WindowEditAnnotation.cpp \
    main.cpp \
    MainWindow.cpp \
    MyCmdLineParser.cpp \
    MyUtils.cpp \
    MyVTKUtils.cpp \
    PanelLayer.cpp \
    PanelPointSet.cpp \
    PanelROI.cpp \
    PanelSurface.cpp \
    PanelVolume.cpp \
    qtcolorpicker.cpp \
    Region2D.cpp \
    Region2DLine.cpp \
    Region2DPolyline.cpp \
    Region2DRectangle.cpp \
    RenderView.cpp \
    RenderView2D.cpp \
    RenderView3D.cpp \
    SurfaceAnnotation.cpp \
    SurfaceLabel.cpp \
    SurfaceOverlay.cpp \
    SurfaceOverlayProperty.cpp \
    SurfaceRegion.cpp \
    SurfaceRegionGroups.cpp \
    TermWidget.cpp \
    ThreadBuildContour.cpp \
    ThreadIOWorker.cpp \
    ToolWindowEdit.cpp \
    ToolWindowMeasure.cpp \
    ToolWindowROIEdit.cpp \
    UIUpdateHelper.cpp \
    VolumeCropper.cpp \
    VolumeFilter.cpp \
    VolumeFilterConvolve.cpp \
    VolumeFilterGradient.cpp \
    VolumeFilterMean.cpp \
    VolumeFilterMedian.cpp \
    VolumeFilterSobel.cpp \
    vtkSimpleLabelEdgeFilter.cpp \
    vtkSimpleLabelEdgeFilter3D.cpp \
    WidgetHistogram.cpp \
    WindowConfigureOverlay.cpp \
    WindowQuickReference.cpp \
    FSTrack.cpp \
    track_io/TrackIO.cpp \
    TrackData.cpp \
    Track.cpp \
    LayerTrack.cpp \
    TrackGroup.cpp \
    PanelTrack.cpp \
    LayerPropertyTrack.cpp \
    DialogSaveVolume.cpp \
    DialogReplaceLabel.cpp \
    LayerVolumeTrack.cpp \
    LayerLandmarks.cpp \
    SurfaceROI.cpp \
    ProgressCallback.cpp \
    MainApplication.cpp \
    DialogRepositionSurface.cpp \
    WindowTimeCourse.cpp \
    WidgetTimeCoursePlot.cpp \
    LayerMRIWorkerThread.cpp \
    DialogLabelStats.cpp \
    VolumeFilterWorkerThread.cpp \
    FSGroupDescriptor.cpp \
    WindowGroupPlot.cpp \
    WidgetGroupPlot.cpp \
    SurfaceSpline.cpp \
    DialogLoadSurfaceOverlay.cpp \
    DialogReloadLayer.cpp \
    DialogSmoothSurface.cpp \
    DialogLineProfile.cpp \
    LayerLineProfile.cpp \
    LayerPropertyLineProfile.cpp \
    LayerConnectomeMatrix.cpp \
    DialogLoadConnectome.cpp \
    PanelConnectomeMatrix.cpp \
    LayerPropertyConnectomeMatrix.cpp \
    PanelAllLayers.cpp \
    DialogLoadSurface.cpp \
    LayerFCD.cpp \
    LayerPropertyFCD.cpp \
    PanelFCD.cpp \
    LayerFCDWorkerThread.cpp \
    DialogLoadFCD.cpp \
    VolumeFilterErode.cpp \
    VolumeFilterDilate.cpp \
    VolumeFilterOpen.cpp \
    VolumeFilterClose.cpp \
    DialogSetCamera.cpp \
    DialogThresholdVolume.cpp \
    DialogVolumeSegmentation.cpp \
    LabelTreeWidget.cpp \
    SplineTreeWidget.cpp \
    DialogLoadTransform.cpp \
    Interactor3DROIEdit.cpp \
    DialogAddPointSetStat.cpp \
    BinaryTreeNode.cpp \
    BinaryTreeEdge.cpp \
    BinaryTreeView.cpp \
    DialogSelectSplines.cpp \
    SurfacePath.cpp \
    Interactor3DPathEdit.cpp \
    DialogCustomFill.cpp \
    DialogSurfaceLabelOperations.cpp \
    geos/GeodesicMatting.cpp \
    geos/kde.cpp \
    GeoSWorker.cpp \
    QVTKWidget/QVTKWidget.cxx \
    QVTKWidget/QVTKPaintEngine.cxx \
    BusyIndicator.cpp \
    vtkInteractorStyleMyTrackballCamera.cxx \
    FlowLayout.cpp \
    WindowLayerInfo.cpp \
    DialogScreenshotOverlay.cpp

HEADERS  += \
    Annotation2D.h \
    BrushProperty.h \
    CommandEdit.h \
    Contour2D.h \
    Cursor2D.h \
    Cursor3D.h \
    CursorFactory.h \
    DialogAbout.h \
    DialogControlPointComment.h \
    DialogCropVolume.h \
    DialogGradientFilter.h \
    DialogLoadDTI.h \
    DialogLoadODF.h \
    DialogLoadPointSet.h \
    DialogLoadVolume.h \
    DialogMovePoint.h \
    DialogNewAnnotation.h \
    DialogPreferences.h \
    DialogNewPointSet.h \
    DialogNewROI.h \
    DialogNewVolume.h \
    DialogSaveAllVolumes.h \
    DialogSavePointSet.h \
    DialogSaveScreenshot.h \
    DialogTransformSurface.h \
    DialogTransformVolume.h \
    DialogVolumeFilter.h \
    DialogWriteMovieFrames.h \
    FloatingStatusBar.h \
    FSLabel.h\
    FSPointSet.h \
    FSSurface.h \
    FSVolume.h \
    GenericRenderView.h \
    InfoTreeWidget.h \
    Interactor.h \
    Interactor2D.h \
    Interactor2DMeasure.h \
    Interactor2DNavigate.h \
    Interactor2DPointSetEdit.h \
    Interactor2DROIEdit.h \
    Interactor2DVolumeCrop.h \
    Interactor2DVolumeEdit.h \
    Interactor2DVoxelEdit.h \
    Interactor3D.h \
    Interactor3DMeasure.h \
    Interactor3DNavigate.h \
    Interactor3DVolumeCrop.h \
    Layer.h \
    LayerCollection.h \
    LayerDTI.h \
    LayerEditable.h \
    LayerODF.h \
    LayerPLabel.h \
    LayerPointSet.h \
    LayerProperty.h \
    LayerPropertyDTI.h \
    LayerPropertyMRI.h \
    LayerPropertyODF.h \
    LayerPropertyPointSet.h \
    LayerPropertyROI.h \
    LayerPropertySurface.h \
    LayerROI.h \
    LayerSurface.h \
    LayerTreeWidget.h \
    LayerVolumeBase.h \
    LayerMRI.h \
    LivewireTool.h \
    LUTDataHolder.h \
    MainWindow.h \
    MigrationDefs.h \
    MyCmdLineParser.h \
    MyUtils.h \
    PanelLayer.h \
    PanelODF.h \
    PanelPointSet.h \
    PanelROI.h \
    PanelSurface.h \
    PanelVolume.h \
    Region3D.h \
    ScribblePromptWorker.h \
    ToolWindowLesionPopup.h \
    TorchScriptModule.h \
    VolumeFilterOptimal.h \
    WindowEditAnnotation.h \
    qtcolorpicker.h \
    Region2D.h \
    Region2DLine.h \
    Region2DPolyline.h \
    Region2DRectangle.h \
    RenderView.h \
    RenderView2D.h \
    RenderView3D.h \
    SurfaceLabel.h \
    SurfaceAnnotation.h \
    SurfaceOverlay.h \
    SurfaceOverlayProperty.h \
    SurfaceRegion.h \
    SurfaceRegionGroups.h \
    TermWidget.h \
    ThreadBuildContour.h \
    ThreadIOWorker.h \
    ToolWindowEdit.h \
    ToolWindowMeasure.h \
    ToolWindowROIEdit.h \
    UIUpdateHelper.h \
    VolumeCropper.h \
    VolumeFilter.h \
    VolumeFilterSobel.h \
    vtkSimpleLabelEdgeFilter.h \
    vtkSimpleLabelEdgeFilter3D.h \
    WidgetHistogram.h \
    WindowConfigureOverlay.h \
    WindowQuickReference.h \
    FSTrack.h \
    TrackData.h \
    LayerTrack.h \
    TrackGroup.h \
    PanelTrack.h \
    LayerPropertyTrack.h \
    DialogSaveVolume.h \
    DialogReplaceLabel.h \
    LayerVolumeTrack.h \
    LayerLandmarks.h \
    SurfaceROI.h \
    ProgressCallback.h \
    MainApplication.h \
    DialogRepositionSurface.h \
    WindowTimeCourse.h \
    WidgetTimeCoursePlot.h \
    LayerMRIWorkerThread.h \
    DialogLabelStats.h \
    VolumeFilterWorkerThread.h \
    FSGroupDescriptor.h \
    WindowGroupPlot.h \
    WidgetGroupPlot.h \
    SurfaceSpline.h \
    DialogLoadSurfaceOverlay.h \
    DialogReloadLayer.h \
    DialogSmoothSurface.h \
    DialogLineProfile.h \
    LayerLineProfile.h \
    LayerPropertyLineProfile.h \
    LayerConnectomeMatrix.h \
    DialogLoadConnectome.h \
    PanelConnectomeMatrix.h \
    LayerPropertyConnectomeMatrix.h \
    PanelAllLayers.h \
    DialogLoadSurface.h \
    LayerFCD.h \
    LayerPropertyFCD.h \
    PanelFCD.h \
    LayerFCDWorkerThread.h \
    DialogLoadFCD.h \
    VolumeFilterErode.h \
    VolumeFilterDilate.h \
    VolumeFilterOpen.h \
    VolumeFilterClose.h \
    DialogSetCamera.h \
    DialogThresholdVolume.h \
    DialogVolumeSegmentation.h \
    LabelTreeWidget.h \
    SplineTreeWidget.h \
    DialogLoadTransform.h \
    Interactor3DROIEdit.h \
    DialogAddPointSetStat.h \
    BinaryTreeNode.h \
    BinaryTreeEdge.h \
    BinaryTreeView.h \
    DialogSelectSplines.h \
    SurfacePath.h \
    Interactor3DPathEdit.h \
    DialogCustomFill.h \
    DialogSurfaceLabelOperations.h \
    geos/GeodesicMatting.h \
    geos/kde.h \
    GeoSWorker.h \
    QVTKWidget/QVTKWidget.h \
    BusyIndicator.h \
    QVTKWidget/QVTKPaintEngine.h \
    vtkInteractorStyleMyTrackballCamera.h \
    FlowLayout.h \
    WindowLayerInfo.h \
    DialogScreenshotOverlay.h

FORMS    += MainWindow.ui \
    DialogControlPointComment.ui \
    DialogLoadODF.ui \
    DialogMovePoint.ui \
    DialogNewAnnotation.ui \
    DialogSaveAllVolumes.ui \
    DialogTransformSurface.ui \
    PanelODF.ui \
    PanelVolume.ui \
    PanelSurface.ui \
    PanelROI.ui \
    PanelPointSet.ui \
    DialogLoadVolume.ui \
    DialogPreferences.ui \
    ToolWindowLesionPopup.ui \
    ToolWindowMeasure.ui \
    ToolWindowEdit.ui \
    DialogNewVolume.ui \
    DialogNewROI.ui \
    DialogNewPointSet.ui \
    ToolWindowROIEdit.ui \
    DialogLoadPointSet.ui \
    DialogLoadDTI.ui \
    WindowConfigureOverlay.ui \
    DialogTransformVolume.ui \
    DialogCropVolume.ui \
    DialogSaveScreenshot.ui \
    DialogVolumeFilter.ui \
    DialogAbout.ui \
    DialogWriteMovieFrames.ui \
    DialogGradientFilter.ui \
    WindowEditAnnotation.ui \
    WindowQuickReference.ui \
    FloatingStatusBar.ui \
    TermWidget.ui \
    PanelTrack.ui \
    DialogSavePointSet.ui \
    DialogSaveVolume.ui \
    DialogReplaceLabel.ui \
    DialogRepositionSurface.ui \
    WindowTimeCourse.ui \
    DialogLabelStats.ui \
    WindowGroupPlot.ui \
    DialogLoadSurfaceOverlay.ui \
    DialogReloadLayer.ui \
    DialogSmoothSurface.ui \
    DialogLineProfile.ui \
    DialogLoadConnectome.ui \
    PanelConnectomeMatrix.ui \
    PanelAllLayers.ui \
    DialogLoadSurface.ui \
    PanelFCD.ui \
    DialogLoadFCD.ui \
    DialogSetCamera.ui \
    DialogThresholdVolume.ui \
    DialogVolumeSegmentation.ui \
    DialogLoadTransform.ui \
    DialogAddPointSetStat.ui \
    DialogSelectSplines.ui \
    DialogCustomFill.ui \
    DialogSurfaceLabelOperations.ui \
    WindowLayerInfo.ui \
    DialogScreenshotOverlay.ui

RESOURCES += \
    freeview.qrc

include ($$PWD/json/qjson.pri)

#LIBS += \
#    -lvtkhdf5_hl -lvtkhdf5 -lLSDyna  -lvtkNetCDF_cxx

QMAKE_CXXFLAGS += -DUNICODE -D_FILE_OFFSET_BITS=64 -D_LARGE_FILES \
                   -DDEVELOPMENT -DHAVE_OPENMP

if(SUPPRESS_WARNINGS) {
  QMAKE_CXXFLAGS_WARN_ON += -Wno-deprecated -Wno-write-strings #-Wno-reorder
}

# set this to your local dev directory
FREESURFER_DEV_DIR = /homes/5/rpwang/freesurfer_dev

# set this to your local install bin directory
# freeview.bin will be copied to that directory
FREESURFER_BIN = /homes/5/rpwang/freesurfer/bin

# for linux
unix {
!macx {
  greaterThan(QT_MAJOR_VERSION, 4): QT += x11extras
  greaterThan(QT_MAJOR_VERSION, 4): QMAKE_CXXFLAGS += -fpermissive

  LIBS += \
    -lvtkverdict -lvtkGraphics -lvtkmetaio -lvtkpng -lvtkzlib \
    -lvtksqlite -lvtkImaging -lvtkFiltering -lvtkCommon -lvtksys \
    -lvtkGenericFiltering -lvtkexoIIc -lvtkNetCDF -lvtkVolumeRendering \
    -lvtkRendering -lvtkftgl -lvtkWidgets -lvtkHybrid -lvtkIO -lvtkDICOMParser -lvtkjpeg \
    -lvtkfreetype -lvtkhdf5 -lvtkhdf5_hl -lvtktiff -lvtkexpat -lLSDyna -lvtkNetCDF_cxx

  INCLUDEPATH += /usr/pubsw/packages/vtk/5.10.1/include/vtk-5.10 \
                 $$FREESURFER_DEV_DIR/include $$FREESURFER_DEV_DIR/vtkutils \
                 /usr/pubsw/packages/mni/current/include \
                 $$FREESURFER_DEV_DIR/lineprof

  ITK_PATH = /usr/pubsw/packages/itk/4.13.0

  QMAKE_CXXFLAGS += -I$$FREESURFER_DEV_DIR/include -I$$ITK_PATH/include/ITK-4.13 \
      -I/usr/pubsw/packages/vxl/current/include/vxl/core \
      -I/usr/pubsw/packages/vxl/current/include/vxl/vcl \
      -I/usr/pubsw/packages/vxl/current/include/vxl/v3p/netlib \
      -I/usr/pubsw/packages/vxl/current/include/vxl/v3p/netlib/opt \
      -I/usr/pubsw/packages/petsc/current/include

  CONFIG(debug, debug|release) {
      QMAKE_CXXFLAGS += -g -O0
  }


  LIBS += -L/usr/pubsw/packages/vtk/5.10.1/lib/vtk-5.10 -L/usr/X11R6/lib \
      -lX11 -lXext -lXt -lSM -lICE -lGLU -ldl \
      -L/usr/pubsw/packages/vxl/current/lib -L$$ITK_PATH/lib/InsightToolkit \
      $$FREESURFER_DEV_DIR/utils/libutils.a $$FREESURFER_DEV_DIR/fsgdf/libfsgdf.a \
      $$FREESURFER_DEV_DIR/vtkutils/libvtkutils.a \
      $$FREESURFER_DEV_DIR/lineprof/liblineprof.a \
      $$FREESURFER_DEV_DIR/hipsstubs/libhipsstubs.a $$FREESURFER_DEV_DIR/vtkutils/libvtkutils.a \
      $$FREESURFER_DEV_DIR/rgb/librgb.a $$FREESURFER_DEV_DIR/unix/libunix.a $$FREESURFER_DEV_DIR/dicom/libdicom.a \
   #   $$FREESURFER_DEV_DIR/jpeg/libjpeg.a $$FREESURFER_DEV_DIR/tiff/libtiff.a $$FREESURFER_DEV_DIR/expat/libexpat.a \
      /usr/pubsw/packages/jpeg/6b/lib/libjpeg.a /usr/pubsw/packages/tiff/3.6.1/lib/libtiff.a /usr/pubsw/packages/expat/2.0.1/lib/libexpat.a \
    $$ITK_PATH/lib/libITKIOSpatialObjects-4.13.a $$ITK_PATH/lib/libITKIOXML-4.13.a $$ITK_PATH/lib/libITKLabelMap-4.13.a \
    $$ITK_PATH/lib/libITKQuadEdgeMesh-4.13.a $$ITK_PATH/lib/libITKOptimizers-4.13.a $$ITK_PATH/lib/libITKPolynomials-4.13.a \
    $$ITK_PATH/lib/libITKBiasCorrection-4.13.a $$ITK_PATH/lib/libITKBioCell-4.13.a $$ITK_PATH/lib/libITKIOBMP-4.13.a \
    $$ITK_PATH/lib/libITKIOBioRad-4.13.a $$ITK_PATH/lib/libITKIOBruker-4.13.a $$ITK_PATH/lib/libITKIOCSV-4.13.a \
    $$ITK_PATH/lib/libITKIOGDCM-4.13.a $$ITK_PATH/lib/libitkgdcmMSFF-4.13.a $$ITK_PATH/lib/libitkgdcmDICT-4.13.a \
    $$ITK_PATH/lib/libitkgdcmIOD-4.13.a $$ITK_PATH/lib/libitkgdcmDSED-4.13.a $$ITK_PATH/lib/libitkgdcmCommon-4.13.a \
  $$ITK_PATH/lib/libitkgdcmjpeg8-4.13.a $$ITK_PATH/lib/libitkgdcmjpeg12-4.13.a $$ITK_PATH/lib/libitkgdcmjpeg16-4.13.a \
    $$ITK_PATH/lib/libitkgdcmopenjp2-4.13.a $$ITK_PATH/lib/libitkgdcmcharls-4.13.a $$ITK_PATH/lib/libitkgdcmuuid-4.13.a \
    $$ITK_PATH/lib/libITKIOGE-4.13.a $$ITK_PATH/lib/libITKIOGIPL-4.13.a $$ITK_PATH/lib/libITKIOHDF5-4.13.a $$ITK_PATH/lib/libITKIOJPEG-4.13.a \
    $$ITK_PATH/lib/libITKIOLSM-4.13.a $$ITK_PATH/lib/libITKIOTIFF-4.13.a $$ITK_PATH/lib/libitktiff-4.13.a $$ITK_PATH/lib/libitkjpeg-4.13.a \
    $$ITK_PATH/lib/libITKIOMINC-4.13.a $$ITK_PATH/lib/libitkminc2-4.13.a $$ITK_PATH/lib/libITKIOMRC-4.13.a $$ITK_PATH/lib/libITKIOMesh-4.13.a \
    $$ITK_PATH/lib/libITKgiftiio-4.13.a $$ITK_PATH/lib/libITKEXPAT-4.13.a $$ITK_PATH/lib/libITKIOMeta-4.13.a $$ITK_PATH/lib/libITKMetaIO-4.13.a \
    $$ITK_PATH/lib/libITKIONIFTI-4.13.a $$ITK_PATH/lib/libITKniftiio-4.13.a $$ITK_PATH/lib/libITKznz-4.13.a $$ITK_PATH/lib/libITKIONRRD-4.13.a \
    $$ITK_PATH/lib/libITKNrrdIO-4.13.a $$ITK_PATH/lib/libITKIOPNG-4.13.a $$ITK_PATH/lib/libitkpng-4.13.a $$ITK_PATH/lib/libITKIOSiemens-4.13.a \
    $$ITK_PATH/lib/libITKIOIPL-4.13.a $$ITK_PATH/lib/libITKIOStimulate-4.13.a $$ITK_PATH/lib/libITKIOTransformHDF5-4.13.a $$ITK_PATH/lib/libitkhdf5_cpp.a \
    $$ITK_PATH/lib/libitkhdf5.a $$ITK_PATH/lib/libitkzlib-4.13.a $$ITK_PATH/lib/libITKIOTransformInsightLegacy-4.13.a $$ITK_PATH/lib/libITKIOTransformMatlab-4.13.a \
    $$ITK_PATH/lib/libITKIOTransformBase-4.13.a $$ITK_PATH/lib/libITKTransformFactory-4.13.a $$ITK_PATH/lib/libITKIOVTK-4.13.a \
    $$ITK_PATH/lib/libITKIOImageBase-4.13.a $$ITK_PATH/lib/libITKKLMRegionGrowing-4.13.a $$ITK_PATH/lib/libITKWatersheds-4.13.a \
    $$ITK_PATH/lib/libITKStatistics-4.13.a $$ITK_PATH/lib/libitkNetlibSlatec-4.13.a $$ITK_PATH/lib/libITKSpatialObjects-4.13.a \
    $$ITK_PATH/lib/libITKMesh-4.13.a $$ITK_PATH/lib/libITKTransform-4.13.a $$ITK_PATH/lib/libITKPath-4.13.a $$ITK_PATH/lib/libITKCommon-4.13.a \
    $$ITK_PATH/lib/libitkdouble-conversion-4.13.a $$ITK_PATH/lib/libitksys-4.13.a $$ITK_PATH/lib/libITKVNLInstantiation-4.13.a $$ITK_PATH/lib/libitkvnl_algo-4.13.a \
     $$ITK_PATH/lib/libitkvnl-4.13.a $$ITK_PATH/lib/libitkv3p_netlib-4.13.a $$ITK_PATH/lib/libitknetlib-4.13.a $$ITK_PATH/lib/libitkvcl-4.13.a \
      #/usr/lib64/libuuid.a
      -luuid -lz -lcrypt -ldl -lpthread \
      /usr/pubsw/packages/mni/1.4/lib/libvolume_io.a -L/usr/pubsw/packages/mni/1.4/lib /usr/pubsw/packages/mni/1.4/lib/libminc.a /usr/pubsw/packages/mni/1.4/lib/libnetcdf.a \
      -lvnl_algo -lvnl -lvcl -lnetlib -lv3p_netlib \
      -L/usr/pubsw/packages/petsc/current/lib -lpetscts -lpetscsnes -lpetscksp \
      -lpetscdm -lpetscmat -lpetscvec -lpetscpich -lfmpich \
      /usr/lib64/liblapack.a /usr/lib64/libblas.a -lgfortran -fopenmp

  TARGET = freeview.bin
  DESTDIR = $$FREESURFER_BIN
  }
}

# for mac
macx {

TARGET = FreeView
RC_FILE = resource/icons/freeview.icns

HEADERS  +=  \
    MacHelper.h

OBJECTIVE_SOURCES += \
    MacHelper.mm

greaterThan(QT_MAJOR_VERSION, 4): QT -= x11extras script

QMAKE_CXXFLAGS += -DVCL_CAN_STATIC_CONST_INIT_FLOAT=0

QMAKE_CXXFLAGS -= -DHAVE_OPENMP

QMAKE_MACOSX_DEPLOYMENT_TARGET = 10.9

LIBS += -lvtkHybrid -lvtkVolumeRendering -lvtkRendering -lvtkIO \
          -lvtkGraphics -lvtkGenericFiltering  \
          -lvtkImaging -lvtkDICOMParser -lvtkFiltering -lvtktiff \
          -lvtkCommon -lvtkftgl -lvtkfreetype -lvtkexpat -lvtkjpeg -lvtkpng -lvtksys \
          -lvtkzlib

FREESURFER_DEV_DIR = /Users/rpwang/freesurfer/dev

INCLUDEPATH += /usr/pubsw/packages/vtk/current/include/vtk-5.6 $$FREESURFER_DEV_DIR/include $$FREESURFER_DEV_DIR/vtkutils \
               "/usr/pubsw/packages/mni/current/include" \
               $$FREESURFER_DEV_DIR/lineprof

LIBS += -L/usr/pubsw/packages/vtk/current/lib/vtk-5.6 -framework OpenGL -ldl -lz -framework ApplicationServices \
    -framework CoreServices -framework cocoa -framework IOKit \
    -L/usr/pubsw/packages/vxl/current/lib -L$$ITK_PATH/lib/InsightToolkit \
    $$FREESURFER_DEV_DIR/lib/libutils.a $$FREESURFER_DEV_DIR/lib/libfsgdf.a \
    $$FREESURFER_DEV_DIR/lib/libhipsstubs.a $$FREESURFER_DEV_DIR/lib/libvtkutils.a \
    $$FREESURFER_DEV_DIR/lib/librgb.a $$FREESURFER_DEV_DIR/lib/libunix.a $$FREESURFER_DEV_DIR/lib/libdicom.a \
    $$FREESURFER_DEV_DIR/lib/libjpeg.a $$FREESURFER_DEV_DIR/lib/libtiff.a $$FREESURFER_DEV_DIR/lib/libexpat.a \
    $$FREESURFER_DEV_DIR/lib/liblineprof.a \
    $$ITK_PATH/lib/InsightToolkit/libITKIO.a $$ITK_PATH/lib/InsightToolkit/libITKAlgorithms.a \
    $$ITK_PATH/lib/InsightToolkit/libITKCommon.a $$ITK_PATH/lib/InsightToolkit/libITKMetaIO.a \
    $$ITK_PATH/lib/InsightToolkit/libITKniftiio.a $$ITK_PATH/lib/InsightToolkit/libITKNrrdIO.a \
    $$ITK_PATH/lib/InsightToolkit/libitkpng.a $$ITK_PATH/lib/InsightToolkit/libitksys.a \
    $$ITK_PATH/lib/InsightToolkit/libitktiff.a $$ITK_PATH/lib/InsightToolkit/libitkv3p_netlib.a \
    $$ITK_PATH/lib/InsightToolkit/libitkzlib.a $$ITK_PATH/lib/InsightToolkit/libitkgdcm.a \
    $$ITK_PATH/lib/InsightToolkit/libitkopenjpeg.a $$ITK_PATH/lib/InsightToolkit/libitkjpeg8.a \
    $$ITK_PATH/lib/InsightToolkit/libitkjpeg12.a $$ITK_PATH/lib/InsightToolkit/libitkjpeg16.a \
    $$ITK_PATH/lib/InsightToolkit/libITKDICOMParser.a -lz -ldl -lpthread \
    /usr/pubsw/packages/mni/current/lib/libvolume_io.a -L/usr/pubsw/packages/mni/current/lib \
    /usr/pubsw/packages/mni/current/lib/libminc.a /usr/pubsw/packages/mni/current/lib/libnetcdf.a \
    -lvnl_algo -lvnl -lvcl -lnetlib -lv3p_netlib \
    -L/usr/pubsw/packages/petsc/current/lib -lpetscts -lpetscsnes -lpetscksp \
    -lpetscdm -lpetscmat -lpetscvec -lpetscpich -lpmpich \
    -framework Accelerate /usr/local/gfortran/lib/libgfortran.a

#LIBS -= -L/usr/X11R6/lib -lX11 -lXext -lXt -lSM -lICE -lGLU -lGL

INCLUDEPATH += /usr/local/lib/python3.10/site-packages/torch/include

LIBS += -L/usr/local/lib/python3.10/site-packages/torch/lib -ltorch -lc10 -ltorch_cpu

DESTDIR = ./
}

OTHER_FILES += \
    resource/QuickRef.html \
    Makefile.am

DISTFILES += \
    CMakeLists.txt
