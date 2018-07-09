project(mri_coreg)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_coreg mri_coreg.c)
target_link_libraries(mri_coreg utils)

install(TARGETS mri_coreg DESTINATION bin)
