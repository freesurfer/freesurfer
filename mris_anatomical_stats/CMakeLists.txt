project(mris_anatomical_stats)
include_directories(${mris_anatomical_stats_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_anatomical_stats_SRCS
mris_anatomical_stats.c
)


add_executable(mris_anatomical_stats ${mris_anatomical_stats_SRCS})
target_link_libraries(mris_anatomical_stats ${FS_LIBS})
install(TARGETS mris_anatomical_stats DESTINATION bin)	

