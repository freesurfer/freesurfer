project(mri_make_density_map)
include_directories(${mri_make_density_map_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_make_density_map_SRCS
mri_make_density_map.c
)


add_executable(mri_make_density_map ${mri_make_density_map_SRCS})
target_link_libraries(mri_make_density_map ${FS_LIBS})
install(TARGETS mri_make_density_map DESTINATION bin)	


