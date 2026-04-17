project(dissection_photo_2)

if(BUILD_GUIS)

  include_directories(${CMAKE_CURRENT_SOURCE_DIR})

  set(SOURCES
    ../WidgetImageView.cpp
    ../cnpy.cpp
    ../exif.cpp
    MaskProcessor.cpp
    ../MyUtils.cpp
    MainWindow.cpp
    ProgressWindow.cpp
    DialogWelcome.cpp
    main.cpp
    QThumbnailWidget.cpp
    MakeThumbnailThread.cpp
    dissection_photo.qrc
  )

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
