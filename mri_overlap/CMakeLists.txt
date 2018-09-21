project(mri_overlap)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_overlap mri_overlap.cpp)
add_help(mri_overlap mri_overlap.help.xml)
target_link_libraries(mri_overlap utilscpp)
install(TARGETS mri_overlap DESTINATION bin) 
