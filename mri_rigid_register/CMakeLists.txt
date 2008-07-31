project(mri_rigid_register)
include_directories(${mri_rigid_register_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_rigid_register_SRCS
mri_rigid_register.c
)


add_executable(mri_rigid_register ${mri_rigid_register_SRCS})
target_link_libraries(mri_rigid_register ${FS_LIBS})
install(TARGETS mri_rigid_register DESTINATION bin)	

