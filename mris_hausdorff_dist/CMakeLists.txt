project(mris_hausdorff_dist)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_hausdorff_dist mris_hausdorff_dist.c)
target_link_libraries(mris_hausdorff_dist utils)

install(TARGETS mris_hausdorff_dist DESTINATION bin)
