project(mri_concatenate_lta)
include_directories(${mri_concatenate_lta_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_concatenate_lta_SRCS
mri_concatenate_lta.c
)


add_executable(mri_concatenate_lta ${mri_concatenate_lta_SRCS})
target_link_libraries(mri_concatenate_lta ${FS_LIBS})
install(TARGETS mri_concatenate_lta DESTINATION bin)	

