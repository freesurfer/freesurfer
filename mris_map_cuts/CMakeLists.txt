project(mris_map_cuts)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_map_cuts mris_map_cuts.c)
target_link_libraries(mris_map_cuts utils)

install(TARGETS mris_map_cuts DESTINATION bin)
