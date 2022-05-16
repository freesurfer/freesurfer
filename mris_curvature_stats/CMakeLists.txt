project(mris_curvature_stats)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_curvature_stats mris_curvature_stats.cpp)
target_link_libraries(mris_curvature_stats utils)

host_os()
if( (HOST_OS MATCHES "CentOS8") OR (HOST_OS MATCHES "Ubuntu20") )
   message(STATUS "mris_curvature_stats test will be skipped on ${HOST_OS}")
else()
   add_test_script(NAME mris_curvature_stats_test SCRIPT test.sh)
endif()

install(TARGETS mris_curvature_stats DESTINATION bin)
