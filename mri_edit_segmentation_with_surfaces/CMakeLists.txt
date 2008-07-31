project(mri_edit_segmentation_with_surfaces)
include_directories(${mri_edit_segmentation_with_surfaces_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_edit_segmentation_with_surfaces_SRCS
mri_edit_segmentation_with_surfaces.c
)


add_executable(mri_edit_segmentation_with_surfaces ${mri_edit_segmentation_with_surfaces_SRCS})
target_link_libraries(mri_edit_segmentation_with_surfaces ${FS_LIBS})
install(TARGETS mri_edit_segmentation_with_surfaces DESTINATION bin)	

