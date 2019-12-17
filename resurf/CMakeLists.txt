project(resurf)

if(ITK_FOUND AND VTK_FOUND)

  include(${ITK_USE_FILE})
  include(${VTK_USE_FILE})

  add_compile_options(-Wno-deprecated -Wno-overloaded-virtual -Wno-inconsistent-missing-override -Wno-unused-variable -Wno-unused-local-typedefs -Wno-reorder)
  include_directories(
    ${FS_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}/Code
    ${CMAKE_SOURCE_DIR}/anatomicuts/Code
    ${CMAKE_SOURCE_DIR}/freeview/track_io
    ${CMAKE_SOURCE_DIR}/freeview
    SYSTEM
    ${ITK_INCLUDE_DIRS}
    ${VTK_INCLUDE_DIRS}
    /space/erebus/2/users/vsiless/code/commong/alglib/src
  )

  set(TRACKIO ../freeview/track_io/TrackIO.cpp)

  add_executable(mris_multimodal mris_multimodal.cxx)
  target_link_libraries(mris_multimodal utils ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  install(TARGETS mris_multimodal DESTINATION bin)
  
  add_executable(mris_profileClustering mris_profileClustering.cxx)
  target_link_libraries(mris_profileClustering utils ${OPENGL_LIBRARIES} ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  install(TARGETS mris_profileClustering DESTINATION bin)
 
  add_executable(mris_surf2vtk mris_surf2vtk.cxx)
  target_link_libraries(mris_surf2vtk utils ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  install(TARGETS mris_surf2vtk DESTINATION bin)
 
  add_executable(mris_image2vtk mris_image2vtk.cxx) # Code/vtkNIFTIImageReader.cxx)
  target_link_libraries(mris_image2vtk utils ${ITK_LIBRARIES} ${VTK_LIBRARIES} )
  install(TARGETS mris_image2vtk DESTINATION bin)
  
  #file(GLOB ALGLIB_SRC /space/erebus/2/users/vsiless/code/commong/alglib/src/*.cpp)
  #add_executable(sigmoidFit sigmoidFit.cxx ${ALGLIB_SRC})
  #target_link_libraries(sigmoidFit utils ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  #install(TARGETS sigmoidFit DESTINATION bin)

  add_executable(mri_vessel_segment mri_vessel_segment.cxx)
  target_link_libraries(mri_vessel_segment utils ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  install(TARGETS mri_vessel_segment DESTINATION bin)

  add_executable(mris_multimodal_surface_placement mris_multimodal_surface_placement.cxx ${TRACKIO})
  target_link_libraries(mris_multimodal_surface_placement utils ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  install(TARGETS mris_multimodal_surface_placement DESTINATION bin)

  add_executable(mris_curvature2image mris_curvature2image.cxx)
  target_link_libraries(mris_curvature2image utils ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  install(TARGETS mris_curvature2image DESTINATION bin)

  add_executable(mris_extract_values mris_extract_values.cxx)
  target_link_libraries(mris_extract_values utils ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  install(TARGETS mris_extract_values DESTINATION bin)

endif()
