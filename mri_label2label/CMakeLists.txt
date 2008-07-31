project(mri_label2label)
include_directories(${mri_label2label_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_label2label_SRCS
mri_label2label.c
)


add_executable(mri_label2label ${mri_label2label_SRCS})
target_link_libraries(mri_label2label ${FS_LIBS})
install(TARGETS mri_label2label DESTINATION bin)	

