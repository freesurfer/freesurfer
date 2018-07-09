project(histo_fix_topology)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(histo_fix_topology histo_fix_topology.c)
target_link_libraries(histo_fix_topology utils)

install(TARGETS histo_fix_topology DESTINATION bin)
