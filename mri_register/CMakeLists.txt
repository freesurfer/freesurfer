project(mri_register)
include_directories(${mri_register_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_register_SRCS
mri_register.c
)


add_executable(mri_register ${mri_register_SRCS})
target_link_libraries(mri_register ${FS_LIBS})
install(TARGETS mri_register DESTINATION bin)	


