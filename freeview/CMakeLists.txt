project(freeview)

message(STATUS "FreeView =======================================")
message(STATUS "FREEVIEW CHECK: BUILD_GUIS=${BUILD_GUIS}")
message(STATUS "FREEVIEW CHECK: VTK_FOUND=${VTK_FOUND}")
message(STATUS "FREEVIEW CHECK: ITK_FOUND=${ITK_FOUND}")

if(TARGET vtkutils)
  message(STATUS "FREEVIEW CHECK: vtkutils target exists")
else()
  message(STATUS "FREEVIEW CHECK: vtkutils target DOES NOT exist")
endif()
message(STATUS "FreeView =======================================")

# if(BUILD_GUIS AND VTK_FOUND AND PETSC_FOUND AND TARGET vtkutils)
if(BUILD_GUIS AND VTK_FOUND AND TARGET vtkutils AND ITK_FOUND)

  include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${FS_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}/vtkutils
    SYSTEM
    ${ITK_INCLUDE_DIRS}
    ${VTK_INCLUDE_DIRS}
  )

  if (PETSC_FOUND AND NOT DISABLE_LINEPROF)
     include_directories(${CMAKE_SOURCE_DIR}/lineprof)
  endif()

  if(HOST_OS MATCHES "Rocky9" OR HOST_OS MATCHES "Ubuntu24")
     include_directories(${VTK_INCLUDE_DIRS})
  else()
     include(${VTK_USE_FILE})
  endif()
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
    Interactor3DPointSetEdit.cpp
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
    LayerODF.cpp
    LayerPropertyODF.cpp
    PanelODF.cpp
    DialogTransformSurface.cpp
    Region3D.cpp
    DialogMovePoint.cpp
    VolumeFilterOptimal.cpp
    DialogSaveAllVolumes.cpp
    ToolWindowLesionPopup.cpp
    DialogLoadODF.cpp
    LayerEditRef.cpp
    DialogGifMaker.cpp
    GifWriterWrapper.cpp
    DialogWayPointsAngle.cpp
    ColorTableItem.cpp
    DialogLabelStatsInRange.cpp
    freeview.qrc
  )

  host_os()

  # petsc distribution may not be available on all platforms
  if(PETSC_FOUND AND NOT DISABLE_LINEPROF)
     set(SOURCES ${SOURCES} DialogLineProfile.cpp LayerLineProfile.cpp)
  endif()

  if(BUILD_SCRIBBLE AND (HOST_OS MATCHES "Rocky8" OR HOST_OS MATCHES "Rocky9"))
     message(STATUS "For HOST_OS=${HOST_OS} Build of freeview with ScribblePrompt enabled")
     set(SOURCES ${SOURCES} TorchScriptModule.cpp ScribblePromptWorker.cpp)
     set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_RPATH}:/space/freesurfer/packages/libtorch/2.3.0.dev20240124+cpu/lib")
  else()
     set(BUILD_SCRIBBLE OFF)
     message(STATUS "For HOST_OS=${HOST_OS} Build of freeview with ScribblePrompt disabled")
  endif()

  if(HOST_OS MATCHES "Ubuntu24")
    # do not strip RPATH in install pass since system java libs may not otherwise be found and use RUNPATH not DT_RPATH
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--disable-new-dtags")
  endif()

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
  elseif(${VTK_MAJOR_VERSION} LESS 9)
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/QVTK8)
    set(SOURCES ${SOURCES}
      QVTK8/vtkEventQtSlotConnect.cxx
      QVTK8/vtkQtConnection.cxx
      QVTK8/QVTKInteractor.cxx
      QVTK8/QVTKInteractorAdapter.cxx
      QVTK8/QVTKWidget.cxx
      QVTK8/QVTKPaintEngine.cxx
      QVTK8/QVTKOpenGLNativeWidget.cxx
      QVTK8/QVTKOpenGLWidget.cxx
      QVTK8/QVTKOpenGLWindow.cxx
    )
  else()
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/QVTK9)
    set(SOURCES ${SOURCES}
      QVTK9/vtkEventQtSlotConnect.cxx
      QVTK9/vtkQtConnection.cxx
      QVTK9/QVTKInteractor.cxx
      QVTK9/QVTKInteractorAdapter.cxx
      QVTK9/QVTKOpenGLNativeWidget.cxx
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

  # this option did not work on darwin_x86_64 (Intel mac)
  if(NOT APPLE)
     add_compile_options(-nostartfiles)
  elseif(APPLE_ARM64)
     add_compile_options(-nostartfiles)
  endif()

  if (DISABLE_LINEPROF)
    # add_definitions(-DDISABLE_LINEPROF)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDISABLE_LINEPROF")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DDISABLE_LINEPROF")
  endif()

  if(Qt6_DIR)
    set(QT_LIBRARIES Qt6::Core Qt6::Widgets)
  elseif(Qt5_DIR)
    set(QT_LIBRARIES Qt5::Core Qt5::Widgets)
  endif()

  if(APPLE)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -framework Cocoa -framework IOKit -framework Accelerate")
  else()
    set(QT_LIBRARIES ${QT_LIBRARIES} Qt5::X11Extras)
  endif()

  add_executable(freeview ${SOURCES})

  if(HOST_OS MATCHES "Ubuntu24" AND UB24_SYSLIBS)
     if(ITK_SYSLIBS)
        set(ITK_LIBRARIES "")
     endif()
     if(VTK_SYSLIBS)
        set(VTK_LIBRARIES "")
     endif()
  endif()

  # petsc and/or profiling may not work/compile on some platforms
  if(PETSC_FOUND AND NOT DISABLE_LINEPROF)
     target_link_libraries(freeview
       nifti
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
       ${QUADMATH_LIBRARIES}
     )
  else()
     target_link_libraries(freeview
       nifti
       vtkutils
       ${VTK_LIBRARIES}
       ${ITK_LIBRARIES}
       utils
       ${QT_LIBRARIES}
       ${OPENGL_LIBRARIES}
       ${LAPACK_LIBRARIES}
       ${BLAS_LIBRARIES}
       ${GFORTRAN_LIBRARIES}
       ${QUADMATH_LIBRARIES}
     )
  endif()

  if(NOT APPLE)
    target_link_libraries(freeview ${GLUT_LIBRARIES} ${X11_LIBRARIES})
  endif()

  if(HOST_OS MATCHES "Ubuntu24")
     if(NOT UBUNTU_ARM64)
        set(UB24_ITK_BASE "/usr/lib/${UBUNTU_GCC_LIBDIR}")
        set(UB24_VTK_BASE "/usr/lib/${UBUNTU_GCC_LIBDIR}")
        set(EXT "5.3.so")
     else()
        # no itk system libs in ubuntu24  arm64
        set(UB24_ITK_BASE "/usr/pubsw/packages/itk/5.4_ubuntu24/lib")
        set(UB24_VTK_BASE "/usr/lib/${UBUNTU_GCC_LIBDIR}")
        set(EXT "5.4.so")
     endif()

     set(UB24_ITK_LIBS libITKCommon-${EXT} libITKgiftiio-${EXT} libITKniftiio-${EXT} libitkv3p_netlib-${EXT} libitkvnl-${EXT} libitkvnl_algo-${EXT})

     set(UB24_VTK_LIBS libvtkCommonComputationalGeometry-9.1.so.1 libvtkCommonCore-9.1.so.1 libvtkCommonDataModel-9.1.so.1 libvtkCommonExecutionModel-9.1.so.1 libvtkCommonMath-9.1.so.1 libvtkCommonMisc-9.1.so.1 libvtkCommonTransforms-9.1.so.1 libvtkFiltersCore-9.1.so.1 libvtkFiltersExtraction-9.1.so.1 libvtkFiltersGeneral-9.1.so.1 libvtkFiltersHybrid-9.1.so.1 libvtkFiltersModeling-9.1.so.1 libvtkFiltersSources-9.1.so.1 libvtkGUISupportQt-9.1.so.1 libvtkIOCore-9.1.so.1 libvtkIOExport-9.1.so.1 libvtkIOGeometry-9.1.so.1 libvtkIOImage-9.1.so.1 libvtkIOLegacy-9.1.so.1 libvtkImagingCore-9.1.so.1 libvtkImagingGeneral-9.1.so.1 libvtkImagingMath-9.1.so.1 libvtkImagingMorphological-9.1.so.1 libvtkInteractionStyle-9.1.so.1 libvtkRenderingAnnotation-9.1.so.1 libvtkRenderingCore-9.1.so.1 libvtkRenderingOpenGL2-9.1.so.1)

     foreach(LIB ${UB24_ITK_LIBS})
        target_link_libraries(freeview ${UB24_ITK_BASE}/${LIB})
     endforeach()
     foreach(LIB ${UB24_VTK_LIBS})
        target_link_libraries(freeview ${UB24_VTK_BASE}/${LIB})
     endforeach()

  endif()

  if(LINK_SHARED_AS_NEEDED)
     target_link_options(freeview PRIVATE "LINKER:--as-needed")
  endif()

  if(APPLE)
    mac_deploy_qt(TARGET freeview BUNDLE Freeview.app PLIST Info.plist ICONS resource/icons/freeview.icns)
    install(PROGRAMS mac_wrapper DESTINATION bin RENAME freeview)
  else()
    install(TARGETS freeview DESTINATION bin)
  endif()

  # Use cmake option -DBUILD_SCRIBBLE=ON for in house rocky8 builds (not for release).
  # Compile extra source files above, change build options and link in 3rd party libtorch.
  # There is a related ifdef in ../packages/dcm2niix/CMakeLists.txt
  if(BUILD_SCRIBBLE)
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSCRIBBLE_PROMPT")
     set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DSCRIBBLE_PROMPT")
     if(APPLE)
        message(FATAL_ERROR "Disabled on MacOS; 3rd party homebrew install may not match torch distribution.")
        # include_directories(/usr/local/lib/python3.10/site-packages/torch/include)
        # link_directories(/usr/local/lib/python3.10/site-packages/torch/lib)
        # set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L/usr/local/lib/python3.10/site-packages/torch/lib -ltorch -lc10 -ltorch_cpu")
     else()
        # Rocky8 Linux only from test above
        set(TORCH_DIR "/space/freesurfer/packages/libtorch/2.3.0.dev20240124+cpu")
        include_directories(${TORCH_DIR}/include)
        link_directories(${TORCH_DIR}/lib)
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${TORCH_DIR}/lib -ltorch -lc10 -ltorch_cpu")
        set(CMAKE_CXX_STANDARD 17)
        target_compile_features(freeview PUBLIC cxx_std_17)
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0")
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++1z")

        # add_custom_command(TARGET freeview POST_BUILD COMMAND
        #   ./custom_rpath.sh ${CMAKE_CURRENT_SOURCE_DIR} ${Qt5_DIR} ${VTK_DIR} ${TORCH_DIR} ${CMAKE_INSTALL_PREFIX}
        #   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        #   COMMENT "Set RPATH and install as ${CMAKE_INSTALL_PREFIX}/freeview --> freeview_scribble_prompt."
        # )
     endif()
  endif()

endif()

