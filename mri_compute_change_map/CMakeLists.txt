project(mri_compute_change_map)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_compute_change_map mri_compute_change_map.cpp)
target_link_libraries(mri_compute_change_map utils)

install(TARGETS mri_compute_change_map DESTINATION bin)
