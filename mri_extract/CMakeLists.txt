project(mri_extract)
include_directories(${mri_extract_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_extract_SRCS
mri_extract.c
)


add_executable(mri_extract ${mri_extract_SRCS})
target_link_libraries(mri_extract ${FS_LIBS})
install(TARGETS mri_extract DESTINATION bin)	

