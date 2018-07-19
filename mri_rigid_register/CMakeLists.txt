project(mri_rigid_register)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_rigid_register mri_rigid_register.c)
target_link_libraries(mri_rigid_register utils)

install(TARGETS mri_rigid_register DESTINATION bin)
