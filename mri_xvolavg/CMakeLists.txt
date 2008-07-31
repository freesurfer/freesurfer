project(mri_xvolavg)
include_directories(${mri_xvolavg_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_xvolavg_SRCS
mri_xvolavg.c
)


add_executable(mri_xvolavg ${mri_xvolavg_SRCS})
target_link_libraries(mri_xvolavg ${FS_LIBS})
install(TARGETS mri_xvolavg DESTINATION bin)	
