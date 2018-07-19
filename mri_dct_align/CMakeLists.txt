project(mri_dct_align)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_dct_align mri_dct_align.c)
target_link_libraries(mri_dct_align utils)
install(TARGETS mri_dct_align DESTINATION bin)

add_executable(mri_dct_align_binary mri_dct_align_binary.c)
target_link_libraries(mri_dct_align_binary utils)
install(TARGETS mri_dct_align_binary DESTINATION bin)
