project(retrospective_correction)

if(BUILD_GUIS AND Qt5_FOUND)

  include_directories(${CMAKE_CURRENT_SOURCE_DIR})

  set(SOURCES
    ../WidgetImageView.cpp
    ../cnpy.cpp
    DialogSelectFolder.cpp
    DialogWelcome.cpp
    MainWindow.cpp
    main.cpp
    retrospective_correction.qrc
  )

  add_executable(retrospective_correction ${SOURCES})

  target_link_libraries(retrospective_correction ${QT_LIBRARIES})

  install(TARGETS retrospective_correction DESTINATION bin)

endif()
