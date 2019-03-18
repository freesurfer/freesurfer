project(mri_fslmat_to_lta)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_fslmat_to_lta mri_fslmat_to_lta.cpp)
target_link_libraries(mri_fslmat_to_lta utils)

install(TARGETS mri_fslmat_to_lta DESTINATION bin)
