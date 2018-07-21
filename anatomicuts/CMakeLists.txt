
include_directories(
	Code
	${CMAKE_SOURCE_DIR}/include 
	${CMAKE_SOURCE_DIR}/freeview/track_io 
	${VTK_INCLUDE_DIR}
	#../trc
)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated") # "-Wmaybe-uninitialized")
add_executable(AnatomiCuts
	../freeview/track_io/TrackIO.cpp
	AnatomiCuts.cxx
)

target_link_libraries(AnatomiCuts
		utils
	${ITK_LIBRARIES}
	${VTK_LIBRARIES}
)


