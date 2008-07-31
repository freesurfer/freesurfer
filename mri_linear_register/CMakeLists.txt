project(mri_linear_register)
include_directories(${mri_linear_register_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_linear_register_SRCS
mri_linear_register.c
)


add_executable(mri_linear_register ${mri_linear_register_SRCS})
target_link_libraries(mri_linear_register ${FS_LIBS})
install(TARGETS mri_linear_register DESTINATION bin)	

