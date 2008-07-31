project(mri_surf2vol)
include_directories(${mri_surf2vol_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_surf2vol_SRCS
mri_surf2vol.c
)


add_executable(mri_surf2vol ${mri_surf2vol_SRCS})
target_link_libraries(mri_surf2vol ${FS_LIBS})
install(TARGETS mri_surf2vol DESTINATION bin)	

