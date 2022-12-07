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

  add_executable(fiducials_calibration ${SOURCES})

  target_link_libraries(fiducials_calibration ${QT_LIBRARIES})

  if(APPLE)
    mac_deploy_qt(TARGET fiducials_calibration BUNDLE fiducials_calibration.app PLIST Info.plist)
    install(PROGRAMS mac_wrapper DESTINATION bin RENAME fiducials_calibration)
  else()
    install(TARGETS fiducials_calibration DESTINATION bin)
  endif()

endif()
