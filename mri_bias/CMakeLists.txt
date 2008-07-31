project(mri_bias)
include_directories(${mri_bias_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_bias_SRCS
mri_apply_bias.c
mri_compute_bias.c
)


add_executable(mri_bias ${mri_bias_SRCS})
target_link_libraries(mri_bias ${FS_LIBS})
install(TARGETS mri_bias DESTINATION bin)	

