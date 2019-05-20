project(mris_distance_transform)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_distance_transform mris_distance_transform.cpp)
target_link_libraries(mris_distance_transform utils)

install(TARGETS mris_distance_transform DESTINATION bin)
