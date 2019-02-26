project(mri_stats2seg)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_stats2seg mri_stats2seg.c)
target_link_libraries(mri_stats2seg utils)
install(TARGETS mri_stats2seg DESTINATION bin)
