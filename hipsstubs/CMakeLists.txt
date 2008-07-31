project(hipsstubs)
include_directories(${hipsstubs_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(hipsstubs_SRCS
hipsrepl.c
hipsstubs.c
)


add_library(hipsstubs STATIC ${hipsstubs_SRCS})

