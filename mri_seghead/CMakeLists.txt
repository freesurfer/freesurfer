project(mri_seghead)
include_directories(${mri_seghead_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_seghead_SRCS
mri_seghead.c
)


add_executable(mri_seghead ${mri_seghead_SRCS})
target_link_libraries(mri_seghead ${FS_LIBS})
install(TARGETS mri_seghead DESTINATION bin)	
