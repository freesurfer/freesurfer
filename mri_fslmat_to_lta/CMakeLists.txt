project(mri_fslmat_to_lta)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_fslmat_to_lta mri_fslmat_to_lta.c)
target_link_libraries(mri_fslmat_to_lta utils)

install(TARGETS mri_fslmat_to_lta DESTINATION bin)
