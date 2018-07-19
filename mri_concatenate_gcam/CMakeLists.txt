project(mri_concatenate_gcam)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_concatenate_gcam mri_concatenate_gcam.c)
target_link_libraries(mri_concatenate_gcam utils)

install(TARGETS mri_concatenate_gcam DESTINATION bin)
