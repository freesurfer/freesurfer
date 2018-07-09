project(mri_vol2vol)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_vol2vol mri_vol2vol.c)
target_link_libraries(mri_vol2vol utils)

install(TARGETS mri_vol2vol DESTINATION bin)
