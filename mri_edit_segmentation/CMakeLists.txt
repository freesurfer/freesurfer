project(mri_edit_segmentation)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_edit_segmentation mri_edit_segmentation.c)
target_link_libraries(mri_edit_segmentation utils)

install(TARGETS mri_edit_segmentation DESTINATION bin)
