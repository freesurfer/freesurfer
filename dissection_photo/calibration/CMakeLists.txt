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

  if(APPLE)
    mac_deploy_qt(TARGET calibration_with_fiducials BUNDLE calibration_with_fiducials.app PLIST Info.plist)
    install(PROGRAMS mac_wrapper DESTINATION bin RENAME calibration_with_fiducials)
  else()
    install(TARGETS calibration_with_fiducials DESTINATION bin)
  endif()

endif()
