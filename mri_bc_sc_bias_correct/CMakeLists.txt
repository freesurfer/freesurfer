project(mri_bc_sc_bias_correct)
include_directories(${mri_bc_sc_bias_correct_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_bc_sc_bias_correct_SRCS
mri_bc_sc_bias_correct.c
)


add_executable(mri_bc_sc_bias_correct ${mri_bc_sc_bias_correct_SRCS})
target_link_libraries(mri_bc_sc_bias_correct ${FS_LIBS})
install(TARGETS mri_bc_sc_bias_correct DESTINATION bin)	

