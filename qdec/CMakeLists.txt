project(qdec)

if(BUILD_QDEC
  AND BUILD_GUIS
  AND OPENGL_FOUND
  AND TCLTKTIXBLT_FOUND
  AND VTK_FOUND
  AND KWWidgets_FOUND)
  
  if(APPLE AND NOT BUILD_OSX_TK_GUIS)
  
    install_symlinks(${CMAKE_SOURCE_DIR}/distribution/osx32b/qdec TYPE programs DESTINATION bin)

  else()

    include_directories(
      ${CMAKE_CURRENT_SOURCE_DIR}
      ${FS_INCLUDE_DIRS}
      ${CMAKE_SOURCE_DIR}/fsgdf
      ${CMAKE_SOURCE_DIR}/vtkfsio
      ${CMAKE_SOURCE_DIR}/vtkutils
      ${CMAKE_SOURCE_DIR}/qdecproject
      ${TCLTKTIXBLT_INCLUDE_DIR}
      ${X11_INCLUDE_DIR}
      ${VTK_INCLUDE_DIRS}
      ${KWWidgets_INCLUDE_DIR}
    )

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated")

    # use vtkWrapTcl to create the tcl-wrapped code for these objects
    vtk_wrap_tcl(vtkKWQdecApp.h vtkKWQdecAppTcl.cxx)
    vtk_wrap_tcl(vtkKWQdecWindow.h vtkKWQdecWindowTcl.cxx)
    vtk_wrap_tcl(vtkKWQdecView.h vtkKWQdecViewTcl.cxx)
    vtk_wrap_tcl_init(QdecLibInit.data QdecLibInit.cxx)

    set(SOURCES
      QdecMain.cxx
      vtkKWQdecApp.cxx
      vtkKWQdecAppTcl.cxx
      vtkKWQdecWindow.cxx
      vtkKWQdecWindowTcl.cxx
      vtkKWQdecView.cxx
      vtkKWQdecViewTcl.cxx
      QdecUtilities.cxx
      QdecLibInit.cxx
    )

    if(APPLE)
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -framework Cocoa -framework Carbon -framework IOKit")
    endif()

    add_executable(qdec ${SOURCES})
    target_link_libraries(qdec
      utils
      fsgdf
      vtkfsio
      vtkutils
      qdecproject
      ${OPENGL_LIBRARIES}
      ${TCLTKTIXBLT_LIBRARIES}
      ${X11_LIBRARIES}
      ${VTK_LIBRARIES}
      ${VTK_TCL_LIBRARIES}
      ${KWWidgets_LIBRARIES}
    )

    install(TARGETS qdec DESTINATION bin)

  endif()
  
  install_symlinks(qdec_splash.png TYPE files DESTINATION lib/images)
  install_symlinks(QdecIcons.txt TYPE files DESTINATION lib/resource)

endif()

install(PROGRAMS run_mris_preproc DESTINATION bin)
