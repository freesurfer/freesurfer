project(mri_hires_register)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_hires_register mri_hires_register.c)
target_link_libraries(mri_hires_register utils)
install(TARGETS mri_hires_register DESTINATION bin)

add_executable(mri_nl_align_binary mri_nl_align_binary.c)
target_link_libraries(mri_nl_align_binary utils)
install(TARGETS mri_nl_align_binary DESTINATION bin)

add_executable(mri_linear_align_binary mri_linear_align_binary.c)
target_link_libraries(mri_linear_align_binary utils)
install(TARGETS mri_linear_align_binary DESTINATION bin)

add_executable(mri_compute_distances mri_compute_distances.c)
target_link_libraries(mri_compute_distances utils)
install(TARGETS mri_compute_distances DESTINATION bin)

add_executable(mri_linear_align mri_linear_align.c)
target_link_libraries(mri_linear_align utils)
install(TARGETS mri_linear_align DESTINATION bin)

add_executable(mri_nl_align mri_nl_align.c)
target_link_libraries(mri_nl_align utils)
install(TARGETS mri_nl_align DESTINATION bin)
