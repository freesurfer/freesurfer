project(mri_vol2vol)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_vol2vol mri_vol2vol.cpp)
target_link_libraries(mri_vol2vol utils)

install(TARGETS mri_vol2vol DESTINATION bin)
