project(calibration)

if(BUILD_GUIS AND Qt5_FOUND)

  include_directories(${CMAKE_CURRENT_SOURCE_DIR})

  set(SOURCES
    ../WidgetImageView.cpp
    ../cnpy.cpp
    DialogWelcome.cpp
    MainWindow.cpp
    main.cpp
    calibration.qrc
  )

  add_executable(calibration_with_fiducials ${SOURCES})

  target_link_libraries(calibration_with_fiducials ${QT_LIBRARIES})

  install(TARGETS calibration_with_fiducials DESTINATION bin)

endif()
