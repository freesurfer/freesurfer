project(mri_aseg_edit_reclassify)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_aseg_edit_reclassify mri_aseg_edit_reclassify.c)
target_link_libraries(mri_aseg_edit_reclassify utils)

install(TARGETS mri_aseg_edit_reclassify DESTINATION bin)
