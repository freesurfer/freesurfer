project(mri_z2p)
include_directories(${mri_z2p_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_z2p_SRCS
mri_z2p.c
)


add_executable(mri_z2p ${mri_z2p_SRCS})
target_link_libraries(mri_z2p ${FS_LIBS})
install(TARGETS mri_z2p DESTINATION bin)	


