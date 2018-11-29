project(mri_refine_seg)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_refine_seg mri_refine_seg.cpp)
add_help(mri_refine_seg mri_refine_seg.help.xml)
target_link_libraries(mri_refine_seg utils utilscpp)
install(TARGETS mri_refine_seg DESTINATION bin) 
