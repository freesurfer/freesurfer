project(mri_segment)
include_directories(${mri_segment_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_segment_SRCS
mri_segment.c
)


add_executable(mri_segment ${mri_segment_SRCS})
target_link_libraries(mri_segment ${FS_LIBS})
install(TARGETS mri_segment DESTINATION bin)	
