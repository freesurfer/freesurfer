project(tkregister2)

# todo: move all tkregister2_cmdl stuff to this directory
add_subdirectory(cmdl)

if(OPENGL_FOUND AND TCLTKTIXBLT_FOUND)

  include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/unix
    ${TCLTKTIXBLT_INCLUDE_DIR}
    ${GLUT_INCLUDE_DIR}/GL
  )

  add_executable(tkregister2 tkregister2.c)
  target_link_libraries(tkregister2 utils ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} ${TCLTKTIXBLT_LIBRARIES})

  set(WRAPCODE "\
#!/bin/tcsh -ef\n\
source $FREESURFER_HOME/bin/tcl_setup\n\
exec $FREESURFER_HOME/tktools/tkregister2.bin $argv"
  )
  install_wrapped(TARGETS tkregister2 DESTINATION tktools WRAPCODE ${WRAPCODE})
  install(PROGRAMS ${CMAKE_INSTALL_PREFIX}/bin/tkregister2 DESTINATION bin RENAME tkregister)
  install(PROGRAMS ${CMAKE_INSTALL_PREFIX}/bin/tkregister2 DESTINATION tktools)
  install(PROGRAMS ${CMAKE_INSTALL_PREFIX}/bin/tkregister2 DESTINATION tktools RENAME tkregister)

endif()
