project(mri_mi)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_mi mri_mi.c)
target_link_libraries(mri_mi utils)

install(TARGETS mri_mi DESTINATION bin)
