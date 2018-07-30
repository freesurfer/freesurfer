project(tkregister2)

# todo: move all tkregister2_cmdl stuff up to this directory
add_subdirectory(cmdl)

if(BUILD_GUIS AND OPENGL_FOUND AND TCLTKTIXBLT_FOUND)

  if(APPLE AND NOT BUILD_OSX_TK_GUIS)
  
    install_symlinks(${CMAKE_SOURCE_DIR}/distribution/Mac_32b/tkregister2.bin TYPE programs DESTINATION bin)

  else()

    include_directories(
      ${FS_INCLUDE_DIRS}
      ${CMAKE_SOURCE_DIR}/unix
      ${TCLTKTIXBLT_INCLUDE_DIR}
      ${GLUT_INCLUDE_DIR}/GL
      ${X11_INCLUDE_DIR}
    )

    add_executable(tkregister2.bin tkregister2.c)
    target_link_libraries(tkregister2.bin
      utils
      ${OPENGL_LIBRARIES}
      ${GLUT_LIBRARIES}
      ${X11_LIBRARIES}
      ${TCLTKTIXBLT_LIBRARIES}
    )

    install(TARGETS  tkregister2.bin DESTINATION tktools)

  endif()

  install(PROGRAMS tkregister2_wrapper DESTINATION bin RENAME tkregister2)
  install(PROGRAMS tkregister2_wrapper DESTINATION tktools RENAME tkregister2)
  install(PROGRAMS tkregister2_wrapper DESTINATION tktools RENAME tkregister)

endif()
