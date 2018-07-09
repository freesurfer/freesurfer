project(mris_surface_stats)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_surface_stats mris_surface_stats.c)
target_link_libraries(mris_surface_stats utils)

install(TARGETS mris_surface_stats DESTINATION bin)
