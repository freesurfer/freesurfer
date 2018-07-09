project(mri_transform)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_transform mri_transform.c)
target_link_libraries(mri_transform utils)

install(TARGETS mri_transform DESTINATION bin)
