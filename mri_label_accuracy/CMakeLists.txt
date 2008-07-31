project(mri_label_accuracy)
include_directories(${mri_label_accuracy_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_label_accuracy_SRCS
mri_label_accuracy.c
)


add_executable(mri_label_accuracy ${mri_label_accuracy_SRCS})
target_link_libraries(mri_label_accuracy ${FS_LIBS})
install(TARGETS mri_label_accuracy DESTINATION bin)	

