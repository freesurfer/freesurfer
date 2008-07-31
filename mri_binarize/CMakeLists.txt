project(mri_binarize)
include_directories(${mri_binarize_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_binarize_SRCS
mri_binarize.c
)


add_executable(mri_binarize ${mri_binarize_SRCS})
target_link_libraries(mri_binarize ${FS_LIBS})
install(TARGETS mri_binarize DESTINATION bin)	
