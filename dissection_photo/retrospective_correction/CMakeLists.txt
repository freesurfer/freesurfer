project(retrospective_correction)

if(BUILD_GUIS)

  include_directories(${CMAKE_CURRENT_SOURCE_DIR})

  set(SOURCES
    ../WidgetImageView.cpp
    ../cnpy.cpp
    ../exif.cpp
    ../MyUtils.cpp
    DialogSelectFolder.cpp
    DialogWelcome.cpp
    MainWindow.cpp
    main.cpp
    retrospective_correction.qrc
  )

  add_executable(retrospective_correction ${SOURCES})

  target_link_libraries(retrospective_correction ${QT_LIBRARIES})

  if(APPLE)
    mac_deploy_qt(TARGET retrospective_correction BUNDLE bin/retrospective_correction.app PLIST Info.plist)
    install(PROGRAMS mac_wrapper DESTINATION bin RENAME retrospective_correction)
  else()
    install(TARGETS retrospective_correction DESTINATION bin)
  endif()

endif()
