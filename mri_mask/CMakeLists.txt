project(mri_mask)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_mask mri_mask.c)
add_help(mri_mask mri_mask.help.xml)
target_link_libraries(mri_mask utils)

install(TARGETS mri_mask DESTINATION bin)
