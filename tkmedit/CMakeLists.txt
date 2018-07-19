project(tkmedit)

if(OPENGL_FOUND AND TCLTKTIXBLT_FOUND)

  include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${TCLTKTIXBLT_INCLUDE_DIR}
    ${TIFF_INCLUDE_DIR}
    ${GLUT_INCLUDE_DIR}/GL
  )

  set(SOURCES
    tkmDisplayArea.c
    tkmedit.c
    tkmFunctionalVolume.c
    tkmMeditWindow.c
  )

  add_executable(tkmedit ${SOURCES})
  target_link_libraries(tkmedit utils fsgdf ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} ${TCLTKTIXBLT_LIBRARIES})

  set(WRAPCODE "\
#!/bin/tcsh -ef\n\
source $FREESURFER_HOME/bin/tcl_setup\n\
exec $FREESURFER_HOME/tktools/tkmedit.bin $argv"
  )
  install_wrapped(TARGETS tkmedit DESTINATION tktools WRAPCODE ${WRAPCODE})
  install(PROGRAMS ${CMAKE_INSTALL_PREFIX}/bin/tkmedit DESTINATION tktools)

endif()
