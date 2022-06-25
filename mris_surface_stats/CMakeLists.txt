project(mris_surface_stats)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_surface_stats mris_surface_stats.cpp)
target_link_libraries(mris_surface_stats utils)

add_test_script(NAME mris_surface_stats_test SCRIPT test.sh)

install(TARGETS mris_surface_stats DESTINATION bin)
