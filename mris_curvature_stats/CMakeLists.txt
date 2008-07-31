project(mris_curvature_stats)
include_directories(${mris_curvature_stats_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_curvature_stats_SRCS
mris_curvature_stats.c
)


add_executable(mris_curvature_stats ${mris_curvature_stats_SRCS})
target_link_libraries(mris_curvature_stats ${FS_LIBS})
install(TARGETS mris_curvature_stats DESTINATION bin)	

