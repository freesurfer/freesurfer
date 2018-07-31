project(mris_longitudinal_surfaces)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_longitudinal_surfaces mris_longitudinal_surfaces.c)
target_link_libraries(mris_longitudinal_surfaces utils)

install(TARGETS mris_longitudinal_surfaces DESTINATION bin)
