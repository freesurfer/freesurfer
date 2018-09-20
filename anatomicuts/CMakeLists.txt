
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
  add_executable(AnatomiCuts AnatomiCuts.cxx ${TRACKIO})
  target_link_libraries(AnatomiCuts utils ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  install(TARGETS AnatomiCuts DESTINATION bin)

#AnatomiCuts
  add_executable(AnatomiCuts_correspondences AnatomiCuts_correspondences.cxx ${TRACKIO})
  target_link_libraries(AnatomiCuts_correspondences utils ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  install(TARGETS AnatomiCuts_correspondences DESTINATION bin)

#TRK tools for header
  add_executable(trk_tools trk_tools.cxx ${TRACKIO})
  target_link_libraries(trk_tools utils ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  install(TARGETS trk_tools DESTINATION bin)

#streamline filtering
  add_executable(streamlineFilter StreamlineFilter.cxx ${TRACKIO})
  target_link_libraries(streamlineFilter utils ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  install(TARGETS streamlineFilter DESTINATION bin)
  
endif()
