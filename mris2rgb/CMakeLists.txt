project(mris2rgb)

if(OPENGL_FOUND AND GLUT_FOUND)

  include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/unix
    ${GLUT_INCLUDE_DIR}/GL
    ${TIFF_INCLUDE_DIR}
  )

  add_executable(mris2rgb mris2rgb.c)
  target_link_libraries(mris2rgb utils unix X11 ${OPENGL_LIBRARIES} ${GLUT_LIBRARY})
  install(TARGETS mris2rgb DESTINATION bin)

endif()