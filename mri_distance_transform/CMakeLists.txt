project(mri_distance_transform)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_distance_transform mri_distance_transform.cpp)
target_link_libraries(mri_distance_transform utils)

install(TARGETS mri_distance_transform DESTINATION bin)
