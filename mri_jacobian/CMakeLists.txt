project(mri_jacobian)
include_directories(${mri_jacobian_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_jacobian_SRCS
mri_jacobian.c
)


add_executable(mri_jacobian ${mri_jacobian_SRCS})
target_link_libraries(mri_jacobian ${FS_LIBS})
install(TARGETS mri_jacobian DESTINATION bin)	

