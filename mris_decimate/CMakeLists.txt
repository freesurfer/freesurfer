if(GTS_FOUND AND GLIB2_FOUND)
  project(mris_decimate)

  include_directories(${CMAKE_SOURCE_DIR}/include ${GTS_INCLUDE_DIR} ${GLIB2_INCLUDE_DIRS})

  add_executable(mris_decimate mris_decimate.cpp)
  add_help(mris_decimate mris_decimate.help.xml)
  target_link_libraries(mris_decimate utils ${GTS_LIBRARY} glib-2.0)

  install(TARGETS mris_decimate DESTINATION bin)
endif()
