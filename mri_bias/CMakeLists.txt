project(mri_bias)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_apply_bias mri_apply_bias.c)
target_link_libraries(mri_apply_bias utils)
install(TARGETS mri_apply_bias DESTINATION bin)

add_executable(mri_compute_bias mri_compute_bias.c)
target_link_libraries(mri_compute_bias utils)
install(TARGETS mri_compute_bias DESTINATION bin)
