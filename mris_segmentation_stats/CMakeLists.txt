project(mris_segmentation_stats)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_segmentation_stats mris_segmentation_stats.c)
target_link_libraries(mris_segmentation_stats utils)

install(TARGETS mris_segmentation_stats DESTINATION bin)
