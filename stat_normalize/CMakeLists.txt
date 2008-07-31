project(stat_normalize)
include_directories(${stat_normalize_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(stat_normalize_SRCS
stat_normalize.c
)


add_executable(stat_normalize ${stat_normalize_SRCS})
target_link_libraries(stat_normalize ${FS_LIBS})
install(TARGETS stat_normalize DESTINATION bin)	

