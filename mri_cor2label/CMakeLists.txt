project(mri_cor2label)
include_directories(${mri_cor2label_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_cor2label_SRCS
mri_cor2label.c
)


add_executable(mri_cor2label ${mri_cor2label_SRCS})
target_link_libraries(mri_cor2label ${FS_LIBS})
install(TARGETS mri_cor2label DESTINATION bin)	

