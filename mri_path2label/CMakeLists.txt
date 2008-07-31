project(mri_path2label)
include_directories(${mri_path2label_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_path2label_SRCS
mri_path2label.c
)


add_executable(mri_path2label ${mri_path2label_SRCS})
target_link_libraries(mri_path2label ${FS_LIBS})
install(TARGETS mri_path2label DESTINATION bin)	

