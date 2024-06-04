project(dissection_photo)

if(BUILD_GUIS)

  include_directories(${CMAKE_CURRENT_SOURCE_DIR})

  set(SOURCES
    ../WidgetImageView.cpp
    ../cnpy.cpp
    ../exif.cpp
    ../connected_components/MaskProcessor.cpp
    ../MyUtils.cpp
    MainWindow.cpp
    ProgressWindow.cpp
    main.cpp
    dissection_photo.qrc
  )

  add_executable(dissection_photo ${SOURCES})

  target_link_libraries(dissection_photo ${QT_LIBRARIES})

  if(APPLE)
    mac_deploy_qt(TARGET dissection_photo BUNDLE bin/dissection_photo PLIST Info.plist)
    install(PROGRAMS mac_wrapper DESTINATION bin RENAME dissection_photo)
  else()
    install(TARGETS dissection_photo DESTINATION bin)
  endif()

endif()
