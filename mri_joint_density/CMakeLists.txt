project(mri_joint_density)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_joint_density mri_joint_density.c)
target_link_libraries(mri_joint_density utils)

install(TARGETS mri_joint_density DESTINATION bin)
