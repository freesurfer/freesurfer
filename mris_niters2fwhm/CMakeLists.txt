project(mris_niters2fwhm)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_niters2fwhm mris_niters2fwhm.c)
target_link_libraries(mris_niters2fwhm utils)

install(TARGETS mris_niters2fwhm DESTINATION bin)
