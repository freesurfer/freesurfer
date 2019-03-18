project(mri_surf_smooth)
cmake_minimum_required(VERSION 3.9)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_surf_smooth mri_surf_smooth.c)
target_link_libraries(mri_surf_smooth utils)

install(TARGETS mri_surf_smooth DESTINATION bin)
