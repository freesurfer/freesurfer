project(mri_fit_bias)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_fit_bias mri_fit_bias.cpp)
target_link_libraries(mri_fit_bias utils)

install(TARGETS mri_fit_bias DESTINATION bin)
