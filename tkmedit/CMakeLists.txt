project(tkmedit)

if(BUILD_GUIS AND OPENGL_FOUND AND TCLTKTIXBLT_FOUND)

  if(APPLE AND NOT BUILD_OSX_TK_GUIS)

    install_symlinks(${CMAKE_SOURCE_DIR}/distribution/osx32b/tkmedit TYPE programs DESTINATION bin)

  else()

    include_directories(
      ${FS_INCLUDE_DIRS}
      ${TCLTKTIXBLT_INCLUDE_DIR}
      ${CMAKE_SOURCE_DIR}/packages/tiff
      ${CMAKE_SOURCE_DIR}/packages/glut
      ${X11_INCLUDE_DIR}
    )

    set(SOURCES
      tkmDisplayArea.c
      tkmedit.c
      tkmFunctionalVolume.c
      tkmMeditWindow.c
    )

    add_executable(tkmedit ${SOURCES})
    target_link_libraries(tkmedit
      utils
      fsgdf
      ${OPENGL_LIBRARIES}
      glut
      ${X11_LIBRARIES}
      ${TCLTKTIXBLT_LIBRARIES}
    )

    install(TARGETS tkmedit DESTINATION bin)

  endif()

endif()
