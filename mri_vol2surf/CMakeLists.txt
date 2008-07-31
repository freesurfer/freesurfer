project(mri_vol2surf)
include_directories(${mri_vol2surf_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_vol2surf_SRCS
mri_vol2surf.c
)


add_executable(mri_vol2surf ${mri_vol2surf_SRCS})
target_link_libraries(mri_vol2surf ${FS_LIBS})
install(TARGETS mri_vol2surf DESTINATION bin)	

