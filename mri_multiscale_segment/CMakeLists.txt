project(mri_multiscale_segment)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_multiscale_segment mri_multiscale_segment.cpp)
target_link_libraries(mri_multiscale_segment utils)

install(TARGETS mri_multiscale_segment DESTINATION bin)
