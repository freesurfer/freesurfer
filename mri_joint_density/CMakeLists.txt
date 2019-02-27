project(mri_joint_density)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_joint_density mri_joint_density.cpp)
target_link_libraries(mri_joint_density utils)

install(TARGETS mri_joint_density DESTINATION bin)
