project(mri_concat)
include_directories(${mri_concat_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_concat_SRCS
mri_concat.c
)


add_executable(mri_concat ${mri_concat_SRCS})
target_link_libraries(mri_concat ${FS_LIBS})
install(TARGETS mri_concat DESTINATION bin)	

