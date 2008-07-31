project(mri_label2vol)
include_directories(${mri_label2vol_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_label2vol_SRCS
mri_label2vol.c
)


add_executable(mri_label2vol ${mri_label2vol_SRCS})
target_link_libraries(mri_label2vol ${FS_LIBS})
install(TARGETS mri_label2vol DESTINATION bin)	

