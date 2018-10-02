project(mris_mesh_subdivide)

if(VTK_FOUND)

  include_directories(${FS_INCLUDE_DIRS} ${VTK_INCLUDE_DIRS})

  add_executable(mris_mesh_subdivide mris_mesh_subdivide.cxx)
  target_link_libraries(mris_mesh_subdivide utils ${VTK_LIBRARIES})

  install(TARGETS mris_mesh_subdivide DESTINATION bin)

endif()
