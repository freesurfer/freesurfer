project(mri_make_bem_surfaces)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_make_bem_surfaces mri_make_bem_surfaces.cpp)
target_link_libraries(mri_make_bem_surfaces utils)

install(TARGETS mri_make_bem_surfaces DESTINATION bin)
