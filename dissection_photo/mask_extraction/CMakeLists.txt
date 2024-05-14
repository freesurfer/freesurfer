project(mask_extraction)

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
    mask_extraction.qrc
  )

  add_executable(mask_extraction ${SOURCES})

  target_link_libraries(mask_extraction ${QT_LIBRARIES})

  if(APPLE)
    mac_deploy_qt(TARGET mask_extraction BUNDLE bin/mask_extraction.app PLIST Info.plist)
    install(PROGRAMS mac_wrapper DESTINATION bin RENAME mask_extraction)
  else()
    install(TARGETS mask_extraction DESTINATION bin)
  endif()

endif()
