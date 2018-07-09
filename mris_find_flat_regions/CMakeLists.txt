project(mris_find_flat_regions)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_find_flat_regions mris_find_flat_regions.c)
target_link_libraries(mris_find_flat_regions utils)

install(TARGETS mris_find_flat_regions DESTINATION bin)
