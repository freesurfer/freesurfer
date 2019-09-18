project(histo_fix_topology)

include_directories(${FS_INCLUDE_DIRS})

add_executable(histo_fix_topology histo_fix_topology.cpp)
target_link_libraries(histo_fix_topology utils)

install(TARGETS histo_fix_topology DESTINATION bin)
