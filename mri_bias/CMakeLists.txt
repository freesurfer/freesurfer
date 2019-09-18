project(mri_bias)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_apply_bias mri_apply_bias.cpp)
target_link_libraries(mri_apply_bias utils)
install(TARGETS mri_apply_bias DESTINATION bin)

add_executable(mri_compute_bias mri_compute_bias.cpp)
target_link_libraries(mri_compute_bias utils)
install(TARGETS mri_compute_bias DESTINATION bin)
