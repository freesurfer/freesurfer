project(anatomicuts)

if(ITK_FOUND AND VTK_FOUND)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated -Wmaybe-uninitialized")
include( ${ITK_USE_FILE} ) 
add_executable(AnatomiCuts
	../freeview/track_io/TrackIO.cpp
	AnatomiCuts.cxx
)
target_link_libraries(AnatomiCuts
	utils
	${ITK_LIBRARIES}
	${VTK_LIBRARIES}
)
	include_directories(
		Code
		${FS_INCLUDE_DIRS}
		${CMAKE_SOURCE_DIR}/freeview/track_io 
		${VTK_INCLUDE_DIRS}
	)

	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated -Wmaybe-uninitialized")


	add_executable(AnatomiCuts_correspondences
		../freeview/track_io/TrackIO.cpp
		AnatomiCuts_correspondences.cxx
	)
	target_link_libraries(AnatomiCuts_correspondences
		utils
		${ITK_LIBRARIES}
		${VTK_LIBRARIES}
	)

endif()
