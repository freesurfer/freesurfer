project(nmovie_qt)

if(BUILD_GUIS AND Qt5_FOUND)

  set(CMAKE_AUTOUIC ON)
  set(CMAKE_AUTOMOC ON)

  include_directories(${CMAKE_CURRENT_SOURCE_DIR})

  set(SOURCES
    main.cpp
    MainWindow.cpp
    RenderWidget.cpp
  )

  add_executable(nmovie_qt ${SOURCES})

  target_link_libraries(nmovie_qt Qt5::Core Qt5::Widgets)

  install(TARGETS nmovie_qt DESTINATION bin)

endif()
