project(mri_ms_fitparms)
include_directories(${mri_ms_fitparms_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_ms_fitparms_SRCS
mri_ms_fitparms.c
)


add_executable(mri_ms_fitparms ${mri_ms_fitparms_SRCS})
target_link_libraries(mri_ms_fitparms ${FS_LIBS})
install(TARGETS mri_ms_fitparms DESTINATION bin)	

