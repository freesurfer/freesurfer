project(mris_surface_stats)
include_directories(${mris_surface_stats_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_surface_stats_SRCS
mris_surface_stats.c
)


add_executable(mris_surface_stats ${mris_surface_stats_SRCS})
target_link_libraries(mris_surface_stats ${FS_LIBS})
install(TARGETS mris_surface_stats DESTINATION bin)	

