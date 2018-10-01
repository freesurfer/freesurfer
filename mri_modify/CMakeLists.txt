project(mri_modify)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_modify mri_modify.cpp)
target_link_libraries(mri_modify utils)
install(TARGETS mri_modify DESTINATION bin)

add_executable(mri_copy_params mri_copy_params.cpp)
target_link_libraries(mri_copy_params utils)
install(TARGETS mri_copy_params DESTINATION bin)
