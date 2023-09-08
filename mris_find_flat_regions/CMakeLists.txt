project(mris_find_flat_regions)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_find_flat_regions mris_find_flat_regions.cpp)
add_help(mris_find_flat_regions mris_find_flat_regions.help.xml)
target_link_libraries(mris_find_flat_regions utils)

install(TARGETS mris_find_flat_regions DESTINATION bin)
