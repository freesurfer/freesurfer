project(mris_multiscale_stats)
include_directories(${mris_multiscale_stats_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_multiscale_stats_SRCS
mris_multiscale_stats.c
)


add_executable(mris_multiscale_stats ${mris_multiscale_stats_SRCS})
target_link_libraries(mris_multiscale_stats ${FS_LIBS})
install(TARGETS mris_multiscale_stats DESTINATION bin)	


