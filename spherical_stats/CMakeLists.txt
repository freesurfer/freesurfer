project(spherical_stats)
include_directories(${spherical_stats_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(spherical_stats_SRCS
spherical_stats.c
)


add_executable(spherical_stats ${spherical_stats_SRCS})
target_link_libraries(spherical_stats ${FS_LIBS})
install(TARGETS spherical_stats DESTINATION bin)	

