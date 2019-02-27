project(mri_rigid_register)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_rigid_register mri_rigid_register.cpp)
target_link_libraries(mri_rigid_register utils)

install(TARGETS mri_rigid_register DESTINATION bin)
