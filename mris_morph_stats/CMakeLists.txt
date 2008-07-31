project(mris_morph_stats)
include_directories(${mris_morph_stats_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_morph_stats_SRCS
mris_morph_stats.c
)


add_executable(mris_morph_stats ${mris_morph_stats_SRCS})
target_link_libraries(mris_morph_stats ${FS_LIBS})
install(TARGETS mris_morph_stats DESTINATION bin)	

