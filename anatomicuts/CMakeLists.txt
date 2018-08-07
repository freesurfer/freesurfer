project(anatomicuts)

if(ITK_FOUND AND VTK_FOUND)

  include_directories(
    Code
    ${FS_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}/freeview/track_io 
    ${VTK_INCLUDE_DIRS}
  )

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated -Wmaybe-uninitialized")

  add_executable(AnatomiCuts
    AnatomiCuts.cxx
    ../freeview/track_io/TrackIO.cpp
  )
  target_link_libraries(AnatomiCuts utils ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  install(TARGETS AnatomiCuts DESTINATION bin)

  add_executable(AnatomiCuts_correspondences
    AnatomiCuts_correspondences.cxx
    ../freeview/track_io/TrackIO.cpp
  )
  target_link_libraries(AnatomiCuts_correspondences utils ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  install(TARGETS AnatomiCuts_correspondences DESTINATION bin)

endif()
