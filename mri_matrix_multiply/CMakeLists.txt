project(mri_matrix_multiply)
include_directories(${mri_matrix_multiply_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_matrix_multiply_SRCS
mri_matrix_multiply.c
)


add_executable(mri_matrix_multiply ${mri_matrix_multiply_SRCS})
target_link_libraries(mri_matrix_multiply ${FS_LIBS})
install(TARGETS mri_matrix_multiply DESTINATION bin)	

