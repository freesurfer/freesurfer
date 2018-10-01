project(mri_probedicom)

if(OPENGL_FOUND)

  include_directories(
    ${FS_INCLUDE_DIRS}
    ${GLUT_INCLUDE_DIR}/GL
    ${X11_INCLUDE_DIR}
  )

  add_executable(mri_probedicom mri_probedicom.c)
  target_link_libraries(mri_probedicom utils ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} ${X11_LIBRARIES})
  install(TARGETS mri_probedicom DESTINATION bin)

endif()
