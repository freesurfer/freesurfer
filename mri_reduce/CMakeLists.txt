project(mri_reduce)
include_directories(${mri_reduce_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_reduce_SRCS
mri_reduce.c
)


add_executable(mri_reduce ${mri_reduce_SRCS})
target_link_libraries(mri_reduce ${FS_LIBS})
install(TARGETS mri_reduce DESTINATION bin)	

