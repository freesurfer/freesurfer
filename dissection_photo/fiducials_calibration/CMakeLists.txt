project(fiducials_calibration)

if(BUILD_GUIS)

  include_directories(${CMAKE_CURRENT_SOURCE_DIR})

  set(SOURCES
    ../WidgetImageView.cpp
    ../cnpy.cpp
    ../exif.cpp
    DialogWelcome.cpp
    MainWindow.cpp
    main.cpp
    fiducials_calibration.qrc
  )

  add_executable(fiducials_calibration ${SOURCES})

  target_link_libraries(fiducials_calibration ${QT_LIBRARIES})

  if(APPLE)
    mac_deploy_qt(TARGET fiducials_calibration BUNDLE bin/fiducials_calibration.app PLIST Info.plist)
    install(PROGRAMS mac_wrapper DESTINATION bin RENAME fiducials_calibration)
  else()
    install(TARGETS fiducials_calibration DESTINATION bin)
  endif()

endif()
