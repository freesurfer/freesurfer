project(mri_copy_values)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_copy_values mri_copy_values.c)
target_link_libraries(mri_copy_values utils)

install(TARGETS mri_copy_values DESTINATION bin)
