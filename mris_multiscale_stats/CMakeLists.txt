project(mris_multiscale_stats)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_multiscale_stats mris_multiscale_stats.cpp)
target_link_libraries(mris_multiscale_stats utils)

install(TARGETS mris_multiscale_stats DESTINATION bin)
