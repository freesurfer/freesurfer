project(mri_fwhm)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_fwhm mri_fwhm.c)
target_link_libraries(mri_fwhm utils)

install(TARGETS mri_fwhm DESTINATION bin)
