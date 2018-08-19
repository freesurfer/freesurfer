
project(resurf)

if(ITK_FOUND AND VTK_FOUND)

  include(${ITK_USE_FILE})  

  include_directories(
    ${FS_INCLUDE_DIRS}
    ${ITK_INCLUDE_DIRS}
    ${VTK_INCLUDE_DIRS}
    Code
    ../anatomicuts/Code/
  )


  add_executable(mris_surf2resurf.bin
    mris_surf2resurf.cxx 
  )
  target_link_libraries(mris_surf2resurf.bin utils ${ITK_LIBRARIES} ${VTK_LIBRARIES} )
  install(TARGETS mris_surf2resurf.bin DESTINATION bin)
  install(PROGRAMS mris_surf2resurf_wrapper DESTINATION bin RENAME mris_surf2resurf)

  
endif()
