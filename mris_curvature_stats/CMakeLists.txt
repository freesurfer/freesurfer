project(mris_curvature_stats)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_curvature_stats mris_curvature_stats.c)
target_link_libraries(mris_curvature_stats utils)

install(TARGETS mris_curvature_stats DESTINATION bin)
