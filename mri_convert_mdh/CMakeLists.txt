project(mri_convert_mdh)
include_directories(${mri_convert_mdh_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_convert_mdh_SRCS
mri_convert_mdh.c
)


add_executable(mri_convert_mdh ${mri_convert_mdh_SRCS})
target_link_libraries(mri_convert_mdh ${FS_LIBS})
install(TARGETS mri_convert_mdh DESTINATION bin)	
