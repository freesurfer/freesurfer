project(mris_make_map_surfaces)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_make_map_surfaces mris_make_map_surfaces.cpp)
add_help(mris_make_map_surfaces mris_make_map_surfaces.help.xml)
target_link_libraries(mris_make_map_surfaces utils)

install(TARGETS mris_make_map_surfaces DESTINATION bin)
