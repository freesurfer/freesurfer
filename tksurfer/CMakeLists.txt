project(tksurfer)

if(BUILD_GUIS AND OPENGL_FOUND AND TCLTKTIXBLT_FOUND)

  include_directories(
    ${FS_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}/unix
    ${TCLTKTIXBLT_INCLUDE_DIR}
    ${TIFF_INCLUDE_DIR}
    ${GLUT_INCLUDE_DIR}/GL
    ${X11_INCLUDE_DIR}
  )

  add_executable(tksurfer.bin tksurfer.c)
  target_link_libraries(tksurfer.bin
    utils
    fsgdf
    ${OPENGL_LIBRARIES}
    ${GLUT_LIBRARIES}
    ${X11_LIBRARIES}
    ${TCLTKTIXBLT_LIBRARIES}
  )
  
  install(TARGETS  tksurfer.bin DESTINATION tktools)
  install(PROGRAMS tksurfer_wrapper DESTINATION bin RENAME tksurfer)
  install(PROGRAMS tksurfer_wrapper DESTINATION tktools RENAME tksurfer)

endif()
