project(mris_fwhm)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_fwhm mris_fwhm.c)
target_link_libraries(mris_fwhm utils)

install(TARGETS mris_fwhm DESTINATION bin)
