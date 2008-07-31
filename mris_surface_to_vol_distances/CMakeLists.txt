project(mris_surface_to_vol_distances)
include_directories(${mris_surface_to_vol_distances_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_surface_to_vol_distances_SRCS
mris_surface_to_vol_distances.c
)


add_executable(mris_surface_to_vol_distances ${mris_surface_to_vol_distances_SRCS})
target_link_libraries(mris_surface_to_vol_distances ${FS_LIBS})
install(TARGETS mris_surface_to_vol_distances DESTINATION bin)	

