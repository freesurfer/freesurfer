project(spherical_stats)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(spherical_stats spherical_stats.c)
target_link_libraries(spherical_stats utils)

install(TARGETS spherical_stats DESTINATION bin)
