project(resurf)

if(ITK_FOUND AND VTK_FOUND)

  include(${ITK_USE_FILE})

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated -Wno-maybe-uninitialized -Wno-unused-variable -Wno-unused-local-typedefs")

  include_directories(
    ${FS_INCLUDE_DIRS}
    ${ITK_INCLUDE_DIRS}
    ${VTK_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}/anatomicuts/Code
    Code
  )
  
  #set (CMAKE_CXX_FLAGS "-fpermissive")
  add_executable(mris_multimodal mris_multimodal.cxx)
  target_link_libraries(mris_multimodal utils ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  install(TARGETS mris_multimodal DESTINATION bin)
  
  add_executable(mris_profileClustering mris_profileClustering.cxx)
  target_link_libraries(mris_profileClustering utils ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  install(TARGETS mris_profileClustering DESTINATION bin)
 
  add_executable(mris_surf2vtk mris_surf2vtk.cxx)
  target_link_libraries(mris_surf2vtk utils ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  install(TARGETS mris_surf2vtk DESTINATION bin)
 
  #add_executable(mris_image2vtk mris_image2vtk.cxx)
  #target_link_libraries(mris_image2vtk utils ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  #install(TARGETS mris_image2vtk DESTINATION bin)
  
endif()
