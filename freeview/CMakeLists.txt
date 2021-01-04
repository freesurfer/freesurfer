project(freeview)

if(BUILD_GUIS AND VTK_FOUND AND PETSC_FOUND AND TARGET vtkutils)

  include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${FS_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}/vtkutils
    ${CMAKE_SOURCE_DIR}/lineprof
    SYSTEM
    ${ITK_INCLUDE_DIRS}
    ${VTK_INCLUDE_DIRS}
  )

  include(${VTK_USE_FILE})
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
    DialogCustomFill.cpp
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
    VolumeFilterBoundary.cpp
    vtkSimpleLabelEdgeFilter.cpp
    vtkSimpleLabelEdgeFilter3D.cpp
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
    DialogSurfaceLabelOperations.cpp
    GeoSWorker.cpp
    geos/GeodesicMatting.cpp
    geos/kde.cpp
    BusyIndicator.cpp
    vtkInteractorStyleMyTrackballCamera.cxx
    FlowLayout.cpp
    WindowLayerInfo.cpp
    WindowEditAnnotation.cpp
    DialogNewAnnotation.cpp
    DialogScreenshotOverlay.cpp
    DialogControlPointComment.cpp
    freeview.qrc
  )

  # build the QVTKWidget utilities
  if(${VTK_MAJOR_VERSION} LESS 7)
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/QVTK5)
    set(SOURCES ${SOURCES} QVTK5/QVTKWidget.cxx QVTK5/QVTKPaintEngine.cxx)
  elseif(${VTK_MAJOR_VERSION} LESS 8)
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/QVTK7)
    set(SOURCES ${SOURCES}
      QVTK7/vtkEventQtSlotConnect.cxx
      QVTK7/vtkQtConnection.cxx
      QVTK7/QVTKInteractor.cxx
      QVTK7/QVTKInteractorAdapter.cxx
      QVTK7/QVTKWidget.cxx
      QVTK7/QVTKPaintEngine.cxx
    )
  endif()

  # fix for retina screens interacting with QVTKWidget
  if(${VTK_MAJOR_VERSION} LESS 8 AND APPLE)
    set(SOURCES ${SOURCES} MacRetina.mm)
  endif()

  # osx helper functions
  if(APPLE)
    set(SOURCES ${SOURCES} MacHelper.mm)
  endif()

  add_compile_options(-Wno-deprecated -Wno-reorder -Wno-uninitialized -Wno-unused-variable -Wno-sign-compare
                      -Wno-overloaded-virtual -Wno-inconsistent-missing-override -Wno-deprecated-declarations)

  if(NOT FREEVIEW_LINEPROF)
    add_definitions(-DDISABLE_LINEPROF)
  endif()

  set(QT_LIBRARIES Qt5::Core Qt5::Widgets)

  if(APPLE)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -framework Cocoa -framework IOKit -framework Accelerate")
  else()
    set(QT_LIBRARIES ${QT_LIBRARIES} Qt5::X11Extras)
  endif()

  add_executable(freeview ${SOURCES})
  target_link_libraries(freeview
    vtkutils
    lineprof
    ${VTK_LIBRARIES}
    ${PETSC_LIBRARIES}
    ${ITK_LIBRARIES}
    utils
    ${QT_LIBRARIES}
    ${OPENGL_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES}
    ${GFORTRAN_LIBRARIES}
  )

  if(NOT APPLE)
    target_link_libraries(freeview glut ${X11_LIBRARIES})
  endif()

  if(APPLE)
    mac_deploy_qt(TARGET freeview BUNDLE Freeview.app PLIST Info.plist ICONS resource/icons/freeview.icns)
    install(PROGRAMS mac_wrapper DESTINATION bin RENAME freeview)
  else()
    install(TARGETS freeview DESTINATION bin)
  endif()

endif()
