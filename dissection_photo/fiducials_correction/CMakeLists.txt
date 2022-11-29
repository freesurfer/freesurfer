project(fiducials_correction)

if(BUILD_GUIS AND Qt5_FOUND)

  include_directories(${CMAKE_CURRENT_SOURCE_DIR})

  set(SOURCES
    ../WidgetImageView.cpp
    ../cnpy.cpp
    MainWindow.cpp
    main.cpp
    fiducials_correction.qrc
  )

  add_executable(fiducials_correction ${SOURCES})

  target_link_libraries(fiducials_correction ${QT_LIBRARIES})

  install(TARGETS fiducials_correction DESTINATION bin)

endif()
