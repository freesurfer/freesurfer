project(mri_fuse_segmentations)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_fuse_segmentations mri_fuse_segmentations.cpp)
add_help(mri_fuse_segmentations mri_fuse_segmentations.help.xml)
target_link_libraries(mri_fuse_segmentations utils)

add_test_script(NAME mri_fuse_segmentations_test SCRIPT test.sh)

install(TARGETS mri_fuse_segmentations DESTINATION bin)
