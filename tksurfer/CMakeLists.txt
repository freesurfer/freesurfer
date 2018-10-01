project(tksurfer)

if(BUILD_GUIS AND OPENGL_FOUND AND TCLTKTIXBLT_FOUND)

  if(APPLE AND NOT BUILD_OSX_TK_GUIS)
  
    install_symlinks(${CMAKE_SOURCE_DIR}/distribution/osx32b/tksurfer TYPE programs DESTINATION bin)

  else()

    include_directories(
      ${FS_INCLUDE_DIRS}
      ${CMAKE_SOURCE_DIR}/unix
      ${TCLTKTIXBLT_INCLUDE_DIR}
      ${TIFF_INCLUDE_DIR}
      ${GLUT_INCLUDE_DIR}/GL
      ${X11_INCLUDE_DIR}
    )

    add_executable(tksurfer tksurfer.c)
    target_link_libraries(tksurfer
      utils
      fsgdf
      ${OPENGL_LIBRARIES}
      ${GLUT_LIBRARIES}
      ${X11_LIBRARIES}
      ${TCLTKTIXBLT_LIBRARIES}
    )
    
    install(TARGETS tksurfer DESTINATION bin)

  endif()

endif()
