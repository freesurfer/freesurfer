project(mris_surface_change)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_surface_change mris_surface_change.c)
target_link_libraries(mris_surface_change utils)

install(TARGETS mris_surface_change DESTINATION bin)
