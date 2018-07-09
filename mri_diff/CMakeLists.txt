project(mri_diff)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_diff mri_diff.c)
target_link_libraries(mri_diff utils)

install(TARGETS mri_diff DESTINATION bin)
