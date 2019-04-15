project(mri_copy_values)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_copy_values mri_copy_values.cpp)
target_link_libraries(mri_copy_values utils)

install(TARGETS mri_copy_values DESTINATION bin)
