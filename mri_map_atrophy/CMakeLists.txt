project(mri_map_atrophy)
include_directories(${mri_map_atrophy_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_map_atrophy_SRCS
mri_map_atrophy.c
)


add_executable(mri_map_atrophy ${mri_map_atrophy_SRCS})
target_link_libraries(mri_map_atrophy ${FS_LIBS})
install(TARGETS mri_map_atrophy DESTINATION bin)	


