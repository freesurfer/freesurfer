project(dissection_photo)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
if(BUILD_GUIS)

  set(CMAKE_AUTOUIC ON)
  set(CMAKE_AUTOMOC ON)
  set(CMAKE_AUTORCC ON)

  if(Qt6_DIR)
    set(QT_LIBRARIES Qt6::Core Qt6::Widgets)
  elseif(Qt5_DIR)
    set(QT_LIBRARIES Qt5::Core Qt5::Widgets)
  endif()

  add_compile_options(-Wno-deprecated -Wno-reorder -Wno-uninitialized -Wno-unused-variable -Wno-sign-compare
                      -Wno-overloaded-virtual -Wno-deprecated-declarations)
                      
                      
  if(APPLE)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -framework Cocoa -framework IOKit -framework Accelerate -lz")
  else()
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--copy-dt-needed-entries")
  endif()

  include_directories(${CMAKE_CURRENT_SOURCE_DIR})
  add_subdirectory(retrospective_correction)
  add_subdirectory(connected_components)
  add_subdirectory(fiducials_calibration)
  add_subdirectory(fiducials_correction)
  add_subdirectory(mask_extraction)
  add_subdirectory(dissection_photo)
  add_subdirectory(nnUNet)

endif()

# should this be in the if statement, is there any reason to install w/o GUI?
install_configured(dissection_photo_gui.sh DESTINATION bin)