project(stim_polar)

if(OPENGL_FOUND)

  include_directories(${FS_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/packages/glut ${X11_INCLUDE_DIR})

  add_executable(stim_polar stim_polar.cpp)
  target_link_libraries(stim_polar utils ${OPENGL_LIBRARIES} glut ${X11_LIBRARIES})
  install(TARGETS stim_polar DESTINATION bin)

endif()
