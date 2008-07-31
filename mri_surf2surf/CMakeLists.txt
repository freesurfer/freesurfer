project(mri_surf2surf)
include_directories(${mri_surf2surf_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_surf2surf_SRCS
mri_surf2surf.c
)


add_executable(mri_surf2surf ${mri_surf2surf_SRCS})
target_link_libraries(mri_surf2surf ${FS_LIBS})
install(TARGETS mri_surf2surf DESTINATION bin)	

