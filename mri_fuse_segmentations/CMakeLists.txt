project(mri_fuse_segmentations)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_fuse_segmentations mri_fuse_segmentations.c)
target_link_libraries(mri_fuse_segmentations utils)

add_test_script(NAME mri_fuse_segmentations_test SCRIPT test.py DEPENDS mri_fuse_segmentations)

install(TARGETS mri_fuse_segmentations DESTINATION bin)
