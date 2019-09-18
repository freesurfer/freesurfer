project(mri_build_priors)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_build_priors mri_build_priors.cpp)
target_link_libraries(mri_build_priors utils)

install(TARGETS mri_build_priors DESTINATION bin)
