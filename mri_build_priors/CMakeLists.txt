project(mri_build_priors)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_build_priors mri_build_priors.c)
target_link_libraries(mri_build_priors utils)

install(TARGETS mri_build_priors DESTINATION bin)
