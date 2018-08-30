
project(anatomicuts)

if(VTK_FOUND)

  include_directories(
    ${FS_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}/Code
    ${CMAKE_SOURCE_DIR}/freeview/track_io
    ${ITK_INCLUDE_DIRS}
    ${VTK_INCLUDE_DIRS}
  )

  include(${ITK_USE_FILE})  

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated -Wmaybe-uninitialized")

  set(TRACKIO ../freeview/track_io/TrackIO.cpp)

#AnatomiCuts correspondences
  add_executable(AnatomiCuts.bin AnatomiCuts.cxx ${TRACKIO})
  target_link_libraries(AnatomiCuts.bin utils ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  install(TARGETS AnatomiCuts.bin DESTINATION bin)
  install(PROGRAMS AnatomiCuts_wrapper DESTINATION bin RENAME AnatomiCuts)

#AnatomiCuts
  add_executable(AnatomiCuts_correspondences.bin AnatomiCuts_correspondences.cxx ${TRACKIO})
  target_link_libraries(AnatomiCuts_correspondences.bin utils ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  install(TARGETS AnatomiCuts_correspondences.bin DESTINATION bin)
  install(PROGRAMS AnatomiCuts_correspondences_wrapper DESTINATION bin RENAME AnatomiCuts_correspondences)

#TRK tools for header
  add_executable(trk_tools.bin trk_tools.cxx ${TRACKIO})
  target_link_libraries(trk_tools.bin utils ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  install(TARGETS trk_tools.bin DESTINATION bin)
  install(PROGRAMS trk_tools_wrapper DESTINATION bin RENAME trk_tools)

#streamline filtering
  add_executable(streamlineFilter.bin StreamlineFilter.cxx ${TRACKIO})
  target_link_libraries(streamlineFilter.bin utils ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  install(TARGETS streamlineFilter.bin DESTINATION bin)
  install(PROGRAMS streamlineFilter_wrapper DESTINATION bin RENAME streamlineFilter)
  
endif()
