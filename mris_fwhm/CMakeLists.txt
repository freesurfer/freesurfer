project(mris_fwhm)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_fwhm mris_fwhm.c)
target_link_libraries(mris_fwhm utils)

install(TARGETS mris_fwhm DESTINATION bin)
