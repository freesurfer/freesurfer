project(mri_edit_wm_with_aseg)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_edit_wm_with_aseg mri_edit_wm_with_aseg.cpp)
add_help(mri_edit_wm_with_aseg mri_edit_wm_with_aseg.help.xml)
target_link_libraries(mri_edit_wm_with_aseg utils)

install(TARGETS mri_edit_wm_with_aseg DESTINATION bin)
