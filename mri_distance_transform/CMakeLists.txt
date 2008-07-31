project(mri_distance_transform)
include_directories(${mri_distance_transform_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_distance_transform_SRCS
mri_distance_transform.cpp
)


add_executable(mri_distance_transform ${mri_distance_transform_SRCS})
target_link_libraries(mri_distance_transform ${FS_LIBS})
install(TARGETS mri_distance_transform DESTINATION bin)	

