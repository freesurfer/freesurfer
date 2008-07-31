project(mri_edit_wm_with_aseg)
include_directories(${mri_edit_wm_with_aseg_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_edit_wm_with_aseg_SRCS
mri_edit_wm_with_aseg.c
)


add_executable(mri_edit_wm_with_aseg ${mri_edit_wm_with_aseg_SRCS})
target_link_libraries(mri_edit_wm_with_aseg ${FS_LIBS})
install(TARGETS mri_edit_wm_with_aseg DESTINATION bin)	

