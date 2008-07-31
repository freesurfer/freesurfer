project(mri_info)
include_directories(${mri_info_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_info_SRCS
mri_info.c
)


add_executable(mri_info ${mri_info_SRCS})
target_link_libraries(mri_info ${FS_LIBS})
install(TARGETS mri_info DESTINATION bin)	

