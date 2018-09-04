project(mris_mesh_subdivide)

if(VTK_FOUND)

  include_directories(${FS_INCLUDE_DIRS} ${VTK_INCLUDE_DIRS})

  add_executable(mris_mesh_subdivide.bin mris_mesh_subdivide.cxx)
  target_link_libraries(mris_mesh_subdivide.bin utils ${VTK_LIBRARIES})

  install(PROGRAMS mris_mesh_subdivide_wrapper DESTINATION bin RENAME mris_mesh_subdivide)
  install(TARGETS mris_mesh_subdivide.bin DESTINATION bin)

endif()
