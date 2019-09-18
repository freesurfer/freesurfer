project(mri_make_density_map)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_make_density_map mri_make_density_map.cpp)
target_link_libraries(mri_make_density_map utils)

install(TARGETS mri_make_density_map DESTINATION bin)
