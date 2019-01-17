project(tksurfer)

if(BUILD_GUIS AND OPENGL_FOUND AND TCLTKTIXBLT_FOUND)

  if(APPLE AND NOT BUILD_OSX_TK_GUIS)
  
    install_symlinks(${CMAKE_SOURCE_DIR}/distribution/osx32b/tksurfer TYPE programs DESTINATION bin)

  else()

    include_directories(
      ${FS_INCLUDE_DIRS}
      ${CMAKE_SOURCE_DIR}/unix
      ${CMAKE_SOURCE_DIR}/packages/tiff
      ${CMAKE_SOURCE_DIR}/packages/glut
      ${TCLTKTIXBLT_INCLUDE_DIR}
      ${X11_INCLUDE_DIR}
    )

    add_executable(tksurfer tksurfer.c)
    target_link_libraries(tksurfer
      utils
      fsgdf
      ${OPENGL_LIBRARIES}
      glut
      ${X11_LIBRARIES}
      ${TCLTKTIXBLT_LIBRARIES}
    )
    
    install(TARGETS tksurfer DESTINATION bin)

  endif()

endif()
