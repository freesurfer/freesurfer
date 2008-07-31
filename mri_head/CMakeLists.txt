project(mri_head)
include_directories(${mri_head_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_head_SRCS
mri_head.c
)


add_executable(mri_head ${mri_head_SRCS})
target_link_libraries(mri_head ${FS_LIBS})
install(TARGETS mri_head DESTINATION bin)	

