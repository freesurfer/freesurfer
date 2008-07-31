project(fslutils)
include_directories(${fslutils_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(fslutils_SRCS
fsl_label2voxel.cpp
)


add_executable(fslutils ${fslutils_SRCS})
target_link_libraries(fslutils ${FS_LIBS})
install(TARGETS fslutils DESTINATION bin)	


