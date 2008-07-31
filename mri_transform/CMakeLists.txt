project(mri_transform)
include_directories(${mri_transform_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_transform_SRCS
mri_transform.c
)


add_executable(mri_transform ${mri_transform_SRCS})
target_link_libraries(mri_transform ${FS_LIBS})
install(TARGETS mri_transform DESTINATION bin)	

