project(mri_extract_label)
include_directories(${mri_extract_label_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_extract_label_SRCS
mri_extract_label.c
)


add_executable(mri_extract_label ${mri_extract_label_SRCS})
target_link_libraries(mri_extract_label ${FS_LIBS})
install(TARGETS mri_extract_label DESTINATION bin)	

