project(mris_curvature_stats)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_curvature_stats mris_curvature_stats.cpp)
target_link_libraries(mris_curvature_stats utils)

add_test_script(NAME mris_curvature_stats_test SCRIPT test.sh DEPENDS mris_curvature_stats)

install(TARGETS mris_curvature_stats DESTINATION bin)
