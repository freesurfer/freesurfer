project(mris_multiscale_stats)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_multiscale_stats mris_multiscale_stats.c)
target_link_libraries(mris_multiscale_stats utils)

install(TARGETS mris_multiscale_stats DESTINATION bin)
