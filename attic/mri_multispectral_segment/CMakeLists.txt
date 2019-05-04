project(mri_multispectral_segment)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_multispectral_segment mri_multispectral_segment.cpp)
target_link_libraries(mri_multispectral_segment utils)

install(TARGETS mri_multispectral_segment DESTINATION bin)
