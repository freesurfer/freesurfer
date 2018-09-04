project(mri_surf2vol)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_surf2vol mri_surf2vol.c)
target_link_libraries(mri_surf2vol utils)

install(TARGETS mri_surf2vol DESTINATION bin)
