
project(anatomicuts)

if(ITK_FOUND AND VTK_FOUND)

  include(${ITK_USE_FILE})  

  include_directories(
    ${ITK_INCLUDE_DIRS}
    ${VTK_INCLUDE_DIRS}
    Code
    ${FS_INCLUDE_DIRS}
    #${FS_INCLUDE_DIRS}/utils/
    ${CMAKE_SOURCE_DIR}/freeview/track_io
  )

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated -Wmaybe-uninitialized")

#AnatomiCuts correspondences
  add_executable(AnatomiCuts_correspondences.bin
    AnatomiCuts_correspondences.cxx ../freeview/track_io/TrackIO.cpp
  )
# ../utils/colortab.c
  target_link_libraries(AnatomiCuts_correspondences.bin ${ITK_LIBRARIES} ${VTK_LIBRARIES} utils)
  install(TARGETS AnatomiCuts_correspondences.bin DESTINATION bin)
  install(PROGRAMS AnatomiCuts_correspondences_wrapper DESTINATION bin RENAME AnatomiCuts_correspondences)


#AnatomiCuts
  add_executable(AnatomiCuts.bin
    AnatomiCuts.cxx ../freeview/track_io/TrackIO.cpp
  )
  target_link_libraries(AnatomiCuts.bin utils ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  install(TARGETS AnatomiCuts.bin DESTINATION bin)
  install(PROGRAMS AnatomiCuts_wrapper DESTINATION bin RENAME AnatomiCuts)


#TRK tools for header
  add_executable(trk_tools.bin
    trk_tools.cxx ../freeview/track_io/TrackIO.cpp
  )
  target_link_libraries(trk_tools.bin utils ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  install(TARGETS trk_tools.bin DESTINATION bin)
  install(PROGRAMS source_vtk_wrapper DESTINATION bin RENAME trk_tools)

  
endif()
