project(mri_compute_change_map)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_compute_change_map mri_compute_change_map.c)
target_link_libraries(mri_compute_change_map utils)

install(TARGETS mri_compute_change_map DESTINATION bin)
