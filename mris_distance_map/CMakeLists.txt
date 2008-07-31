project(mris_distance_map)
include_directories(${mris_distance_map_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_distance_map_SRCS
mris_distance_map.c
)


add_executable(mris_distance_map ${mris_distance_map_SRCS})
target_link_libraries(mris_distance_map ${FS_LIBS})
install(TARGETS mris_distance_map DESTINATION bin)	

