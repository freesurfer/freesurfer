project(mri_dct_align)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_dct_align mri_dct_align.cpp)
target_link_libraries(mri_dct_align utils)
install(TARGETS mri_dct_align DESTINATION bin)

add_executable(mri_dct_align_binary mri_dct_align_binary.cpp)
target_link_libraries(mri_dct_align_binary utils)
install(TARGETS mri_dct_align_binary DESTINATION bin)
