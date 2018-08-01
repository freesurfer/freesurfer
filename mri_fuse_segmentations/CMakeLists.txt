project(mri_fuse_segmentations)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_fuse_segmentations mri_fuse_segmentations.c)
target_link_libraries(mri_fuse_segmentations utils)

install(TARGETS mri_fuse_segmentations DESTINATION bin)
