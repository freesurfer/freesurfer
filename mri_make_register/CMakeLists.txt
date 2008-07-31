project(mri_make_register)
include_directories(${mri_make_register_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_make_register_SRCS
mri_make_register.c
)


add_executable(mri_make_register ${mri_make_register_SRCS})
target_link_libraries(mri_make_register ${FS_LIBS})
install(TARGETS mri_make_register DESTINATION bin)	

