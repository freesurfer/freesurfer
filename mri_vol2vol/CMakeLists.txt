project(mri_vol2vol)
include_directories(${mri_vol2vol_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_vol2vol_SRCS
mri_vol2vol.c
)


add_executable(mri_vol2vol ${mri_vol2vol_SRCS})
target_link_libraries(mri_vol2vol ${FS_LIBS})
install(TARGETS mri_vol2vol DESTINATION bin)	


