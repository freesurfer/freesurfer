project(mri_surfcluster)
include_directories(${mri_surfcluster_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_surfcluster_SRCS
mri_surfcluster.c
)


add_executable(mri_surfcluster ${mri_surfcluster_SRCS})
target_link_libraries(mri_surfcluster ${FS_LIBS})
install(TARGETS mri_surfcluster DESTINATION bin)	

