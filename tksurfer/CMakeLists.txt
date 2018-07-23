project(tksurfer)

if(BUILD_GUIS AND OPENGL_FOUND AND TCLTKTIXBLT_FOUND)

  include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/unix
    ${TCLTKTIXBLT_INCLUDE_DIR}
    ${TIFF_INCLUDE_DIR}
    ${GLUT_INCLUDE_DIR}/GL
  )

  add_executable(tksurfer tksurfer.c)
  target_link_libraries(tksurfer utils fsgdf ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} ${TCLTKTIXBLT_LIBRARIES})
  
  set(WRAPCODE "\
#!/bin/tcsh -ef\n\
source $FREESURFER_HOME/bin/tcl_setup\n\
exec $FREESURFER_HOME/tktools/tksurfer.bin $argv"
  )
  install_wrapped(TARGETS tksurfer DESTINATION tktools WRAPCODE ${WRAPCODE})
  install(PROGRAMS ${CMAKE_INSTALL_PREFIX}/bin/tksurfer DESTINATION tktools)

endif()
