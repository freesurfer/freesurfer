project(mri_label_volume)
include_directories(${mri_label_volume_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_label_volume_SRCS
mri_label_volume.c
)


add_executable(mri_label_volume ${mri_label_volume_SRCS})
target_link_libraries(mri_label_volume ${FS_LIBS})
install(TARGETS mri_label_volume DESTINATION bin)	


