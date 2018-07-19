project(mris_topo_fixer)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_topo_fixer mris_topo_fixer.cpp)
target_link_libraries(mris_topo_fixer utilscpp)
install(TARGETS mris_topo_fixer DESTINATION bin)
