project(mri_em_register)
include_directories(${mri_em_register_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_em_register_SRCS
mri_em_register.c
)


add_executable(mri_em_register ${mri_em_register_SRCS})
target_link_libraries(mri_em_register ${FS_LIBS})
install(TARGETS mri_em_register DESTINATION bin)	

