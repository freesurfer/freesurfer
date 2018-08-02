project(mris2rgb)

if(OPENGL_FOUND)

  include_directories(
    ${FS_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}/unix
    ${GLUT_INCLUDE_DIR}/GL
    ${X11_INCLUDE_DIR}
    ${TIFF_INCLUDE_DIR}
  )

  add_executable(mris2rgb mris2rgb.c)
  target_link_libraries(mris2rgb utils ${GLUT_LIBRARIES} ${X11_LIBRARIES} ${OPENGL_LIBRARIES})
  install(TARGETS mris2rgb DESTINATION bin)

endif()