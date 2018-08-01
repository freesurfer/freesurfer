project(mri_surf2surf)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_surf2surf mri_surf2surf.c)
target_link_libraries(mri_surf2surf utils)

install(TARGETS mri_surf2surf DESTINATION bin)
