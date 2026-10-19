project(dissection_photo_2)

if(BUILD_GUIS)

  include_directories(${CMAKE_CURRENT_SOURCE_DIR})

  set(SOURCES
    ../WidgetImageView.cpp
    ../cnpy.cpp
    ../exif.cpp
    AppMotion.cpp
    ReconPage.cpp
    StepperBar.cpp
    MaskProcessor.cpp
    SegmentedControl.cpp
    CompactSlider.cpp
    ImageGridView.cpp
    ../MyUtils.cpp
    AppTheme.cpp
    MainWindow.cpp
    ProgressWindow.cpp
    DialogWelcome.cpp
    DialogProjectsHub.cpp
    DialogSegmentationSettings.cpp
    ProjectCardWidget.cpp
    ProjectManager.cpp
    main.cpp
    QThumbnailWidget.cpp
    MakeThumbnailThread.cpp
    dissection_photo.qrc
  )

  # App icon in Finder and the Dock, named by CFBundleIconFile in Info.plist
  # (main.cpp also sets it at runtime).
  if(APPLE)
    set(APP_ICON resource/dissection_photo.icns)
    set_source_files_properties(${APP_ICON} PROPERTIES MACOSX_PACKAGE_LOCATION Resources)
    list(APPEND SOURCES ${APP_ICON})
  endif()

  add_executable(dissection_photo_2 ${SOURCES})

  target_link_libraries(dissection_photo_2 ${QT_LIBRARIES})

  if(LINK_SHARED_AS_NEEDED)
     target_link_options(dissection_photo_2 PRIVATE "LINKER:--as-needed")
  endif()

  if(APPLE)
    mac_deploy_qt(TARGET dissection_photo_2 BUNDLE bin/dissection_photo_2 PLIST Info.plist)
    install(PROGRAMS mac_wrapper DESTINATION bin RENAME dissection_photo_gui.sh)
  else()
    install(TARGETS dissection_photo_2 DESTINATION bin)
    install_configured(dissection_photo_gui.sh DESTINATION bin)
  endif()

endif()
