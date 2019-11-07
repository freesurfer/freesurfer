project(mri_dist_surf_label)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_dist_surf_label mri_dist_surf_label.cpp)
target_link_libraries(mri_dist_surf_label utils)

install(TARGETS mri_dist_surf_label DESTINATION bin)
