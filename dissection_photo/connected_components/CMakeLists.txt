project(connected_components)

if(BUILD_GUIS AND Qt5_FOUND)

  include_directories(${CMAKE_CURRENT_SOURCE_DIR})

  set(SOURCES
    ../WidgetImageView.cpp
    ../cnpy.cpp
    DialogSelectFolder.cpp
    DialogWelcome.cpp
    MainWindow.cpp
    MaskProcessor.cpp
    main.cpp
    connected_components.qrc
  )

  add_executable(connected_components ${SOURCES})

  target_link_libraries(connected_components ${QT_LIBRARIES})

  install(TARGETS connected_components DESTINATION bin)

endif()
