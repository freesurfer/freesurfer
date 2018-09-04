project(mri_jacobian)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_jacobian mri_jacobian.c)
target_link_libraries(mri_jacobian utils)

install(TARGETS mri_jacobian DESTINATION bin)
