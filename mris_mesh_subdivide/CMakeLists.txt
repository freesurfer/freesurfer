project(mris_mesh_subdivide)


if(VTK_FOUND)

  include_directories(${CMAKE_SOURCE_DIR}/include ${VTK_INCLUDE_DIR})

  add_executable(mris_mesh_subdivide mris_mesh_subdivide.cxx)
  target_link_libraries(mris_mesh_subdivide utils ${VTK_LIBRARIES})

  set(WRAPCODE "\
#!/bin/tcsh -ef\n\
source $FREESURFER_HOME/bin/vtk_setup\n\
if (-e /bin/pwd) then\n\
  setenv HOME `cd $HOME && /bin/pwd`\n\
  setenv SUBJECTS_DIR `cd $SUBJECTS_DIR && /bin/pwd`\n\
endif\n\
mris_mesh_subdivide.bin $argv"
  )
  install_wrapped(TARGETS mris_mesh_subdivide DESTINATION bin WRAPCODE ${WRAPCODE})

endif()
