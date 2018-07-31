project(mris_surface_stats)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_surface_stats mris_surface_stats.c)
target_link_libraries(mris_surface_stats utils)

install(TARGETS mris_surface_stats DESTINATION bin)
