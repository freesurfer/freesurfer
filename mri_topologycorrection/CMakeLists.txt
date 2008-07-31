project(mri_topologycorrection)
include_directories(${mri_topologycorrection_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_topologycorrection_SRCS
mri_topologycorrection.c
)


add_executable(mri_topologycorrection ${mri_topologycorrection_SRCS})
target_link_libraries(mri_topologycorrection ${FS_LIBS})
install(TARGETS mri_topologycorrection DESTINATION bin)	
