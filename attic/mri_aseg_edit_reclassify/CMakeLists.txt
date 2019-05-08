project(mri_aseg_edit_reclassify)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_aseg_edit_reclassify mri_aseg_edit_reclassify.cpp)
target_link_libraries(mri_aseg_edit_reclassify utils)

install(TARGETS mri_aseg_edit_reclassify DESTINATION bin)
