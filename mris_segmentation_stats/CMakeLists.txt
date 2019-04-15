project(mris_segmentation_stats)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_segmentation_stats mris_segmentation_stats.cpp)
target_link_libraries(mris_segmentation_stats utils)

install(TARGETS mris_segmentation_stats DESTINATION bin)
