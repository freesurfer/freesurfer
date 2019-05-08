project(resurf)

if(ITK_FOUND AND VTK_FOUND)

  include(${ITK_USE_FILE})

  add_compile_options(-Wno-deprecated -Wno-overloaded-virtual -Wno-inconsistent-missing-override -Wno-unused-variable -Wno-unused-local-typedefs -Wno-reorder)
  include_directories(
    ${FS_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}/Code
    ${CMAKE_SOURCE_DIR}/anatomicuts/Code
    ${CMAKE_SOURCE_DIR}/freeview
    SYSTEM
    ${ITK_INCLUDE_DIRS}
    ${VTK_INCLUDE_DIRS}
    /space/erebus/2/users/vsiless/code/commong/alglib/src
 )
 file(GLOB ALGLIB_SRC
    /space/erebus/2/users/vsiless/code/commong/alglib/src/*.cpp
)
  file(GLOB FREEVIEW_SRC
    ${CMAKE_SOURCE_DIR}/freeview/FSVolume.cpp
)


    add_executable(mris_multimodal mris_multimodal.cxx)
  target_link_libraries(mris_multimodal utils ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  install(TARGETS mris_multimodal DESTINATION bin)
  
  add_executable(mris_profileClustering mris_profileClustering.cxx)
  target_link_libraries(mris_profileClustering utils ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  install(TARGETS mris_profileClustering DESTINATION bin)
 
  add_executable(mris_surf2vtk mris_surf2vtk.cxx)
  target_link_libraries(mris_surf2vtk utils ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  install(TARGETS mris_surf2vtk DESTINATION bin)
 
  add_executable(mris_image2vtk mris_image2vtk.cxx  ${FREEVIEW_SRC}) # Code/vtkNIFTIImageReader.cxx)
  target_link_libraries(mris_image2vtk utils ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${QT_LIBRARIES})
  install(TARGETS mris_image2vtk DESTINATION bin)
  
  add_executable(sigmoidFit sigmoidFit.cxx ${ALGLIB_SRC})
  target_link_libraries(sigmoidFit utils ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${QTLibraries})
  install(TARGETS sigmoidFit DESTINATION bin)
endif()
