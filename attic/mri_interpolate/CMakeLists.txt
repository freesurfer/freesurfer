project(mri_interpolate)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_interpolate mri_interpolate.cpp)
target_link_libraries(mri_interpolate utils)

install(TARGETS mri_interpolate DESTINATION bin)
