project(mri_fuse_intensity_images)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_fuse_intensity_images mri_fuse_intensity_images.cpp)
target_link_libraries(mri_fuse_intensity_images utils)

install(TARGETS mri_fuse_intensity_images DESTINATION bin)
