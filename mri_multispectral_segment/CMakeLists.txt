project(mri_multispectral_segment)
include_directories(${mri_multispectral_segment_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_multispectral_segment_SRCS
mri_multispectral_segment.c
)


add_executable(mri_multispectral_segment ${mri_multispectral_segment_SRCS})
target_link_libraries(mri_multispectral_segment ${FS_LIBS})
install(TARGETS mri_multispectral_segment DESTINATION bin)	

