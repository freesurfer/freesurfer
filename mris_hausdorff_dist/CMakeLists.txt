project(mris_hausdorff_dist)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_hausdorff_dist mris_hausdorff_dist.cpp)
target_link_libraries(mris_hausdorff_dist utils)

install(TARGETS mris_hausdorff_dist DESTINATION bin)
