project(mri_transform_to_COR)
include_directories(${mri_transform_to_COR_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_transform_to_COR_SRCS
mri_transform_to_COR.c
)


add_executable(mri_transform_to_COR ${mri_transform_to_COR_SRCS})
target_link_libraries(mri_transform_to_COR ${FS_LIBS})
install(TARGETS mri_transform_to_COR DESTINATION bin)	

