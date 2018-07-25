project(qdec)

project(tksurfer)

if(BUILD_GUIS
  AND OPENGL_FOUND
  AND TCLTKTIXBLT_FOUND
  AND VTK_FOUND
  AND KWWidgets_FOUND)

  include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/fsgdf
    ${CMAKE_SOURCE_DIR}/vtkfsio
    ${CMAKE_SOURCE_DIR}/vtkutils
    ${CMAKE_SOURCE_DIR}/qdecproject
    ${TCLTKTIXBLT_INCLUDE_DIR}
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

  add_executable(qdec ${SOURCES})
  target_link_libraries(qdec
    utils
    fsgdf
    vtkfsio
    vtkutils
    qdecproject
    ${TCLTKTIXBLT_LIBRARIES}
    ${OPENGL_LIBRARIES} X11
    ${VTK_LIBRARIES} ${VTK_TCL_LIBRARIES}
    ${KWWidgets_LIBRARIES}
  )
  
  set(WRAPCODE "\
#!/bin/tcsh -ef\n\
source $FREESURFER_HOME/bin/tcl_setup\n\
source $FREESURFER_HOME/bin/vtk_setup\n\
source $FREESURFER_HOME/bin/kww_setup\n\
if (-e /bin/pwd) then\n\
  setenv HOME `cd $HOME && /bin/pwd`\n\
  setenv SUBJECTS_DIR `cd $SUBJECTS_DIR && /bin/pwd`\n\
endif\n\
qdec.bin $argv"
  )
  install_wrapped(TARGETS qdec DESTINATION bin WRAPCODE ${WRAPCODE})

  install_symlinks(qdec_splash.png TYPE files DESTINATION lib/images)
  install_symlinks(QdecIcons.txt TYPE files DESTINATION lib/resource)

endif()

install(PROGRAMS run_mris_preproc DESTINATION bin)
