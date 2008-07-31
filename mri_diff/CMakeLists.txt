project(mri_diff)
include_directories(${mri_diff_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_diff_SRCS
mri_diff.c
)


add_executable(mri_diff ${mri_diff_SRCS})
target_link_libraries(mri_diff ${FS_LIBS})
install(TARGETS mri_diff DESTINATION bin)	
