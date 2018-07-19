project(mris_longitudinal_surfaces)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_longitudinal_surfaces mris_longitudinal_surfaces.c)
target_link_libraries(mris_longitudinal_surfaces utils)

install(TARGETS mris_longitudinal_surfaces DESTINATION bin)
