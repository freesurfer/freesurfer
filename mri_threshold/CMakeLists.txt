project(mri_threshold)
include_directories(${mri_threshold_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_threshold_SRCS
mri_threshold.c
)


add_executable(mri_threshold ${mri_threshold_SRCS})
target_link_libraries(mri_threshold ${FS_LIBS})
install(TARGETS mri_threshold DESTINATION bin)	

