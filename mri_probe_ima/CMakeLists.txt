project(mri_probe_ima)
include_directories(${mri_probe_ima_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_probe_ima_SRCS
mri_probe_ima.c
)


add_executable(mri_probe_ima ${mri_probe_ima_SRCS})
target_link_libraries(mri_probe_ima ${FS_LIBS})
install(TARGETS mri_probe_ima DESTINATION bin)	
