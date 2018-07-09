project(mri_fuse_intensity_images)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_fuse_intensity_images mri_fuse_intensity_images.c)
target_link_libraries(mri_fuse_intensity_images utils)

install(TARGETS mri_fuse_intensity_images DESTINATION bin)
