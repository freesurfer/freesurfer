project(mri_probedicom)

if(OPENGL_FOUND)

  include_directories(${CMAKE_SOURCE_DIR}/include ${GLUT_INCLUDE_DIR}/GL)

  add_executable(mri_probedicom mri_probedicom.c)
  target_link_libraries(mri_probedicom utils ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})
  install(TARGETS mri_probedicom DESTINATION bin)

endif()
