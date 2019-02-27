project(mris_distance_map)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_distance_map mris_distance_map.cpp)
target_link_libraries(mris_distance_map utils)

install(TARGETS mris_distance_map DESTINATION bin)
