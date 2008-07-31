project(mris_transform)
include_directories(${mris_transform_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_transform_SRCS
mris_transform.c
)


add_executable(mris_transform ${mris_transform_SRCS})
target_link_libraries(mris_transform ${FS_LIBS})
install(TARGETS mris_transform DESTINATION bin)	


