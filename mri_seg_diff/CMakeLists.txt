project(mri_seg_diff)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_seg_diff mri_seg_diff.c)
target_link_libraries(mri_seg_diff utils)

install(TARGETS mri_seg_diff DESTINATION bin)
