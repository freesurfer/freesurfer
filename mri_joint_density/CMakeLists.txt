project(mri_joint_density)
include_directories(${mri_joint_density_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_joint_density_SRCS
mri_joint_density.c
)


add_executable(mri_joint_density ${mri_joint_density_SRCS})
target_link_libraries(mri_joint_density ${FS_LIBS})
install(TARGETS mri_joint_density DESTINATION bin)	

