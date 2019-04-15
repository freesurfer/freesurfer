project(mris_niters2fwhm)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_niters2fwhm mris_niters2fwhm.cpp)
target_link_libraries(mris_niters2fwhm utils)

install(TARGETS mris_niters2fwhm DESTINATION bin)
