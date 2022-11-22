project(dissection_photo)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
if(BUILD_GUIS AND Qt5_FOUND)

  set(CMAKE_AUTOUIC ON)
  set(CMAKE_AUTOMOC ON)
  set(CMAKE_AUTORCC ON)
  set(QT_LIBRARIES Qt5::Core Qt5::Widgets)

  add_compile_options(-Wno-deprecated -Wno-reorder -Wno-uninitialized -Wno-unused-variable -Wno-sign-compare
                      -Wno-overloaded-virtual -Wno-inconsistent-missing-override -Wno-deprecated-declarations)
                      
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--copy-dt-needed-entries")
                      
  if(APPLE)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -framework Cocoa -framework IOKit -framework Accelerate")
  else()
    set(QT_LIBRARIES ${QT_LIBRARIES} Qt5::X11Extras)
  endif()

  include_directories(${CMAKE_CURRENT_SOURCE_DIR})
  add_subdirectory(retrospective_correction)
  add_subdirectory(connected_components)
  add_subdirectory(calibration)
  add_subdirectory(fiducials_correction)

endif()
