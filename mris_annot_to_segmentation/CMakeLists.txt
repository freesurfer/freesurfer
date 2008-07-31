project(mris_annot_to_segmentation)
include_directories(${mris_annot_to_segmentation_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_annot_to_segmentation_SRCS
mris_annot_to_segmentation.c
)


add_executable(mris_annot_to_segmentation ${mris_annot_to_segmentation_SRCS})
target_link_libraries(mris_annot_to_segmentation ${FS_LIBS})
install(TARGETS mris_annot_to_segmentation DESTINATION bin)	

