project(mri_fuse_segmentations)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_fuse_segmentations mri_fuse_segmentations.c)
target_link_libraries(mri_fuse_segmentations utils)

install(TARGETS mri_fuse_segmentations DESTINATION bin)
