project(mri_segment)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_segment mri_segment.c)
add_help(mri_segment mri_segment.help.xml)
target_link_libraries(mri_segment utils)

install(TARGETS mri_segment DESTINATION bin)
