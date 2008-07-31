project(mri_copy_values)
include_directories(${mri_copy_values_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_copy_values_SRCS
mri_copy_values.c
)


add_executable(mri_copy_values ${mri_copy_values_SRCS})
target_link_libraries(mri_copy_values ${FS_LIBS})
install(TARGETS mri_copy_values DESTINATION bin)	

