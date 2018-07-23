project(freeview)

if(BUILD_GUIS AND Qt5_FOUND AND VTK_FOUND AND ITK_FOUND AND PETSC_FOUND AND TARGET vtkutils)

  include_directories(
    ${CMAKE_SOURCE_DIR}/freeview
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/vtkutils
    ${CMAKE_SOURCE_DIR}/lineprof
    ${ITK_INCLUDE_DIRS}
    ${VTK_INCLUDE_DIRS}
    ${QT_INCLUDE_DIRS}
  )

  set(CMAKE_AUTOUIC ON)
  set(CMAKE_AUTOMOC ON)
  set(CMAKE_AUTORCC ON)

  set(SOURCES
    Annotation2D.cpp
    CommandEdit.cpp
    Contour2D.cpp
    Cursor2D.cpp
    Cursor3D.cpp
    CursorFactory.cpp
    BrushProperty.cpp 
    DialogAbout.cpp
    DialogCropVolume.cpp
    DialogGradientFilter.cpp
    DialogLoadDTI.cpp
    DialogLoadPointSet.cpp
    DialogLoadVolume.cpp
    DialogNewROI.cpp
    DialogNewPointSet.cpp
    DialogNewVolume.cpp
    DialogPreferences.cpp
    DialogSavePointSet.cpp
    DialogSaveScreenshot.cpp
    DialogTransformVolume.cpp
    DialogVolumeFilter.cpp
    DialogWriteMovieFrames.cpp
    FloatingStatusBar.cpp
    FSLabel.cpp
    FSPointSet.cpp
    FSSurface.cpp
    FSVolume.cpp
    GenericRenderView.cpp
    InfoTreeWidget.cpp
    Interactor.cpp
    Interactor2D.cpp
    Interactor2DMeasure.cpp
    Interactor2DNavigate.cpp
    Interactor2DPointSetEdit.cpp
    Interactor2DROIEdit.cpp
    Interactor2DVolumeCrop.cpp
    Interactor2DVolumeEdit.cpp
    Interactor2DVoxelEdit.cpp
    Interactor3D.cpp
    Interactor3DMeasure.cpp
    Interactor3DNavigate.cpp
    Interactor3DVolumeCrop.cpp
    Interactor3DROIEdit.cpp
    Interactor3DPathEdit.cpp
    Layer.cpp
    LayerCollection.cpp
    LayerDTI.cpp
    LayerEditable.cpp
    LayerMRI.cpp
    LayerPLabel.cpp
    LayerPointSet.cpp
    LayerProperty.cpp
    LayerPropertyDTI.cpp
    LayerPropertyMRI.cpp
    LayerPropertyPointSet.cpp
    LayerPropertyROI.cpp
    LayerPropertySurface.cpp
    LayerROI.cpp
    LayerSurface.cpp
    LayerTreeWidget.cpp
    LayerVolumeBase.cpp
    LabelTreeWidget.cpp
    SplineTreeWidget.cpp
    LivewireTool.cpp
    LUTDataHolder.cpp
    main.cpp
    MyCmdLineParser.cpp
    MyUtils.cpp
    MyVTKUtils.cpp
    PanelLayer.cpp
    PanelPointSet.cpp
    PanelROI.cpp
    PanelSurface.cpp
    PanelVolume.cpp
    qtcolorpicker.cpp
    QVTKWidget.cxx
    QVTKPaintEngine.cxx
    Region2D.cpp
    Region2DLine.cpp
    Region2DPolyline.cpp
    Region2DRectangle.cpp
    RenderView.cpp
    RenderView2D.cpp
    RenderView3D.cpp
    SurfaceAnnotation.cpp
    SurfaceLabel.cpp
    SurfaceOverlay.cpp
    SurfaceOverlayProperty.cpp
    SurfaceRegion.cpp
    SurfaceRegionGroups.cpp
    TermWidget.cpp
    ThreadBuildContour.cpp
    ThreadIOWorker.cpp
    ToolWindowEdit.cpp
    ToolWindowMeasure.cpp
    ToolWindowROIEdit.cpp
    UIUpdateHelper.cpp
    VolumeCropper.cpp
    VolumeFilter.cpp
    VolumeFilterConvolve.cpp
    VolumeFilterGradient.cpp
    VolumeFilterMean.cpp
    VolumeFilterMedian.cpp
    VolumeFilterSobel.cpp
    VolumeFilterErode.cpp
    VolumeFilterDilate.cpp
    VolumeFilterOpen.cpp
    VolumeFilterClose.cpp
    vtkSimpleLabelEdgeFilter.cpp
    WidgetHistogram.cpp
    WindowConfigureOverlay.cpp
    WindowQuickReference.cpp
    FSTrack.cpp
    track_io/TrackIO.cpp
    TrackData.cpp
    Track.cpp
    LayerTrack.cpp
    TrackGroup.cpp
    PanelTrack.cpp
    LayerPropertyTrack.cpp
    DialogSaveVolume.cpp
    MainWindow.cpp
    DialogReplaceLabel.cpp
    LayerVolumeTrack.cpp
    LayerLandmarks.cpp
    SurfaceROI.cpp
    SurfacePath.cpp
    MainApplication.cpp
    ProgressCallback.cpp
    DialogRepositionSurface.cpp
    WindowTimeCourse.cpp
    WidgetTimeCoursePlot.cpp
    LayerMRIWorkerThread.cpp
    DialogLabelStats.cpp
    VolumeFilterWorkerThread.cpp
    FSGroupDescriptor.cpp
    WindowGroupPlot.cpp
    WidgetGroupPlot.cpp
    SurfaceSpline.cpp
    DialogLoadSurfaceOverlay.cpp
    DialogReloadLayer.cpp
    DialogSmoothSurface.cpp
    LayerLineProfile.cpp
    DialogLineProfile.cpp
    LayerPropertyLineProfile.cpp
    LayerConnectomeMatrix.cpp
    LayerPropertyConnectomeMatrix.cpp
    PanelConnectomeMatrix.cpp
    DialogLoadConnectome.cpp
    PanelAllLayers.cpp
    DialogLoadSurface.cpp
    LayerFCD.cpp
    LayerPropertyFCD.cpp
    PanelFCD.cpp
    DialogLoadFCD.cpp
    LayerFCDWorkerThread.cpp
    DialogSetCamera.cpp
    DialogThresholdVolume.cpp
    DialogVolumeSegmentation.cpp
    DialogLoadTransform.cpp
    DialogAddPointSetStat.cpp
    BinaryTreeEdge.cpp
    BinaryTreeNode.cpp
    BinaryTreeView.cpp
    DialogSelectSplines.cpp
    freeview.qrc
  )

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated -Wno-reorder -Wno-uninitialized -Wno-unused-variable -Wno-sign-compare")

  add_executable(freeview ${SOURCES})
  target_link_libraries(freeview 
    Qt5::Core
    Qt5::Widgets
    Qt5::X11Extras
    ${OPENGL_LIBRARIES}
    ${GLUT_LIBRARIES}
    utils
    fsgdf
    vtkutils
    vtkfsio
    lineprof  
    ${VTK_LIBRARIES}
    ${PETSC_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES}
    ${GFORTRAN_LIBRARIES}
    ${ITK_LIBRARIES}
  )

  set(WRAPCODE "\
#!/bin/tcsh -ef\n\
source $FREESURFER_HOME/bin/tcl_setup\n\
source $FREESURFER_HOME/bin/vtk_setup\n\
source $FREESURFER_HOME/bin/qt_setup\n\
if($#argv > 0) then\n\
  if(\\\"$argv[1]\\\" == \\\"-tkmedit\\\") then\n\
    tkmeditfv $argv[2-$#argv]\n\
    exit 0\n\
  endif\n\
endif\n\
freeview.bin $argv:q"
  )
  install_wrapped(TARGETS freeview DESTINATION bin WRAPCODE ${WRAPCODE})

endif()
