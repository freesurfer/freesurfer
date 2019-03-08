project(mri_hausdorff_dist)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_hausdorff_dist mri_hausdorff_dist.cpp)
target_link_libraries(mri_hausdorff_dist utils)

install(TARGETS mri_hausdorff_dist DESTINATION bin)
