project(mri_hausdorff_dist)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_hausdorff_dist mri_hausdorff_dist.c)
target_link_libraries(mri_hausdorff_dist utils)

install(TARGETS mri_hausdorff_dist DESTINATION bin)
