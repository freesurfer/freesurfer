project(log)
include_directories(${log_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(log_SRCS
congraph.c
connect.c
h_logz.c
lmedian.c
map.c
runfuncs.c
)


add_library(log STATIC ${log_SRCS})
