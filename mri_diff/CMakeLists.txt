project(mri_diff)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_diff mri_diff.c)
target_link_libraries(mri_diff utils)

install(TARGETS mri_diff DESTINATION bin)
