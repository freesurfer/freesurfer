project(mri_seg_overlap)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_seg_overlap mri_seg_overlap.cpp)
add_help(mri_seg_overlap mri_seg_overlap.help.xml)
target_link_libraries(mri_seg_overlap utils)
install(TARGETS mri_seg_overlap DESTINATION bin) 
