project(mri_make_bem_surfaces)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_make_bem_surfaces mri_make_bem_surfaces.c)
target_link_libraries(mri_make_bem_surfaces utils)

install(TARGETS mri_make_bem_surfaces DESTINATION bin)
