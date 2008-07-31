project(mri_annotation2label)
include_directories(${mri_annotation2label_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_annotation2label_SRCS
mri_annotation2label.c
)


add_executable(mri_annotation2label ${mri_annotation2label_SRCS})
target_link_libraries(mri_annotation2label ${FS_LIBS})
install(TARGETS mri_annotation2label DESTINATION bin)	

