project(mri_ca_register)
include_directories(${mri_ca_register_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_ca_register_SRCS
mri_ca_register.c
)


add_executable(mri_ca_register ${mri_ca_register_SRCS})
target_link_libraries(mri_ca_register ${FS_LIBS})
install(TARGETS mri_ca_register DESTINATION bin)	

