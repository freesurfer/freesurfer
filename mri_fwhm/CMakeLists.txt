project(mri_fwhm)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_fwhm mri_fwhm.c)
target_link_libraries(mri_fwhm utils)

install(TARGETS mri_fwhm DESTINATION bin)
