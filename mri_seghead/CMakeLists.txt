project(mri_seghead)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_seghead mri_seghead.c)
target_link_libraries(mri_seghead utils)

install(TARGETS mri_seghead DESTINATION bin)
