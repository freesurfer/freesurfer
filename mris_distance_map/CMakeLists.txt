project(mris_distance_map)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_distance_map mris_distance_map.c)
target_link_libraries(mris_distance_map utils)

install(TARGETS mris_distance_map DESTINATION bin)
