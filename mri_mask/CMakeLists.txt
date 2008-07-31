project(mri_mask)
include_directories(${mri_mask_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_mask_SRCS
mri_mask.c
)


add_executable(mri_mask ${mri_mask_SRCS})
target_link_libraries(mri_mask ${FS_LIBS})
install(TARGETS mri_mask DESTINATION bin)	

