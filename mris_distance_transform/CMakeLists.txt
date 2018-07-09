project(mris_distance_transform)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_distance_transform mris_distance_transform.c)
target_link_libraries(mris_distance_transform utils)

install(TARGETS mris_distance_transform DESTINATION bin)
