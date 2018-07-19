project(mri_compute_seg_overlap)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_compute_seg_overlap mri_compute_seg_overlap.c)
add_help(mri_compute_seg_overlap mri_compute_seg_overlap.help.xml)
target_link_libraries(mri_compute_seg_overlap utils)

install(TARGETS mri_compute_seg_overlap DESTINATION bin)
