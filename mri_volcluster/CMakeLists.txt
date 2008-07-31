project(mri_volcluster)
include_directories(${mri_volcluster_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_volcluster_SRCS
mri_volcluster.c
)


add_executable(mri_volcluster ${mri_volcluster_SRCS})
target_link_libraries(mri_volcluster ${FS_LIBS})
install(TARGETS mri_volcluster DESTINATION bin)	

