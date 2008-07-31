project(mri_ca_label)
include_directories(${mri_ca_label_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_ca_label_SRCS
mri_ca_label.c
)


add_executable(mri_ca_label ${mri_ca_label_SRCS})
target_link_libraries(mri_ca_label ${FS_LIBS})
install(TARGETS mri_ca_label DESTINATION bin)	

