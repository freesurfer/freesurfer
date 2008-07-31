project(mri_convert)
include_directories(${mri_convert_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_convert_SRCS
mri_convert.c
)


add_executable(mri_convert ${mri_convert_SRCS})
target_link_libraries(mri_convert ${FS_LIBS})
install(TARGETS mri_convert DESTINATION bin)	
