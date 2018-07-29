project(tkmedit)

if(BUILD_GUIS AND OPENGL_FOUND AND TCLTKTIXBLT_FOUND)

  if(APPLE)
  
    install_symlinks(${CMAKE_SOURCE_DIR}/distribution/Mac_32b/tkmedit.bin TYPE programs DESTINATION bin)

  else()

    include_directories(
      ${FS_INCLUDE_DIRS}
      ${TCLTKTIXBLT_INCLUDE_DIR}
      ${TIFF_INCLUDE_DIR}
      ${GLUT_INCLUDE_DIR}/GL
      ${X11_INCLUDE_DIR}
    )

    set(SOURCES
      tkmDisplayArea.c
      tkmedit.c
      tkmFunctionalVolume.c
      tkmMeditWindow.c
    )

    add_executable(tkmedit.bin ${SOURCES})
    target_link_libraries(tkmedit.bin
      utils
      fsgdf
      ${OPENGL_LIBRARIES}
      ${GLUT_LIBRARIES}
      ${X11_LIBRARIES}
      ${TCLTKTIXBLT_LIBRARIES}
    )

    install(TARGETS  tkmedit.bin DESTINATION tktools)

  endif()

  install(PROGRAMS tkmedit_wrapper DESTINATION bin RENAME tkmedit)
  install(PROGRAMS tkmedit_wrapper DESTINATION tktools RENAME tkmedit)

endif()
