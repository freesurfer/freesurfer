project(mris2rgb)

if(OPENGL_FOUND)

  include_directories(
    ${FS_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}/packages/tiff
    ${CMAKE_SOURCE_DIR}/packages/glut
    ${X11_INCLUDE_DIR}
  )

  add_executable(mris2rgb mris2rgb.cpp oglutil.cpp texfont.cpp)
  target_link_libraries(mris2rgb utils ${OPENGL_LIBRARIES} glut ${X11_LIBRARIES})
  install(TARGETS mris2rgb DESTINATION bin)

endif()
