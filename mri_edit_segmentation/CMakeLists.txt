project(mri_edit_segmentation)
include_directories(${mri_edit_segmentation_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_edit_segmentation_SRCS
mri_edit_segmentation.c
)


add_executable(mri_edit_segmentation ${mri_edit_segmentation_SRCS})
target_link_libraries(mri_edit_segmentation ${FS_LIBS})
install(TARGETS mri_edit_segmentation DESTINATION bin)	

