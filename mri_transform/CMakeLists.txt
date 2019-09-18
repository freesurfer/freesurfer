project(mri_transform)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_transform mri_transform.cpp)
target_link_libraries(mri_transform utils)

install(TARGETS mri_transform DESTINATION bin)
