project(mri_synthesize)
include_directories(${mri_synthesize_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_synthesize_SRCS
mri_synthesize.c
)


add_executable(mri_synthesize ${mri_synthesize_SRCS})
target_link_libraries(mri_synthesize ${FS_LIBS})
install(TARGETS mri_synthesize DESTINATION bin)	



