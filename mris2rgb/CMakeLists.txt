project(mris2rgb)

if(OPENGL_FOUND)

  include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/unix
    ${GLUT_INCLUDE_DIR}/GL
    ${TIFF_INCLUDE_DIR}
  )

  add_executable(mris2rgb mris2rgb.c)
  target_link_libraries(mris2rgb utils ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})
  install(TARGETS mris2rgb DESTINATION bin)

endif()