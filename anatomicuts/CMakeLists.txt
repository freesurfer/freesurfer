project(anatomicuts)

if(ITK_FOUND AND VTK_FOUND)

  include_directories(
    Code
    ${FS_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}/freeview/track_io
    ${ITK_INCLUDE_DIRS}
    ${VTK_INCLUDE_DIRS}
  )

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated -Wmaybe-uninitialized")

  add_executable(AnatomiCuts.bin
    AnatomiCuts.cxx ../freeview/track_io/TrackIO.cpp
  )
  target_link_libraries(AnatomiCuts.bin utils ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  install(TARGETS AnatomiCuts.bin DESTINATION bin)
  install(PROGRAMS AnatomiCuts_wrapper DESTINATION bin RENAME AnatomiCuts)

  add_executable(AnatomiCuts_correspondences.bin
    AnatomiCuts_correspondences.cxx ../freeview/track_io/TrackIO.cpp
  )
  target_link_libraries(AnatomiCuts_correspondences.bin utils ${ITK_LIBRARIES} ${VTK_LIBRARIES})
  install(TARGETS AnatomiCuts_correspondences.bin DESTINATION bin)
  install(PROGRAMS AnatomiCuts_correspondences_wrapper DESTINATION bin RENAME AnatomiCuts_correspondences)

endif()
