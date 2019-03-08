project(mris_decimate)

if(GTS_FOUND)

  include_directories(${FS_INCLUDE_DIRS} ${GTS_INCLUDE_DIRS})

  if(APPLE)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -framework CoreFoundation")
  endif()

  add_executable(mris_decimate main.cpp mris_decimate.cpp)
  add_help(mris_decimate mris_decimate.help.xml)
  target_link_libraries(mris_decimate utils ${GTS_LIBRARIES})

  install(TARGETS mris_decimate DESTINATION bin)

  add_test_script(NAME mris_decimate_test SCRIPT test.sh DEPENDS mris_decimate)

endif()
