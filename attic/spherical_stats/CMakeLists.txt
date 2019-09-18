project(spherical_stats)

include_directories(${FS_INCLUDE_DIRS})

add_executable(spherical_stats spherical_stats.cpp)
target_link_libraries(spherical_stats utils)

install(TARGETS spherical_stats DESTINATION bin)
