project(mri_interpolate)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_interpolate mri_interpolate.c)
target_link_libraries(mri_interpolate utils)

install(TARGETS mri_interpolate DESTINATION bin)
