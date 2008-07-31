project(mri_label_vals)
include_directories(${mri_label_vals_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_label_vals_SRCS
mri_label_vals.c
)


add_executable(mri_label_vals ${mri_label_vals_SRCS})
target_link_libraries(mri_label_vals ${FS_LIBS})
install(TARGETS mri_label_vals DESTINATION bin)	

