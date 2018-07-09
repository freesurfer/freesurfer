project(mri_concatenate_lta)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_concatenate_lta mri_concatenate_lta.c)
target_link_libraries(mri_concatenate_lta utils)

install(TARGETS mri_concatenate_lta DESTINATION bin)
