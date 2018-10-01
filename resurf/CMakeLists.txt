project(resurf)

if(ITK_FOUND AND VTK_FOUND)

  include(${ITK_USE_FILE})  

  include_directories(
    ${FS_INCLUDE_DIRS}
    ${ITK_INCLUDE_DIRS}
    ${VTK_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}/anatomicuts/Code
    Code
  )

  add_executable(mris_surf2resurf mris_surf2resurf.cxx)
  target_link_libraries(mris_surf2resurf utils ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  install(TARGETS mris_surf2resurf DESTINATION bin)
  
endif()
