project(qdecproject)

include_directories(${qdecproject_SOURCE_DIR}
${qdecproject_SOURCE_DIR}/include
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

FILE (GLOB sfiles *.c* )

SET(qdecproject_SRCS
${sfiles}
)

add_library(qdecproject STATIC ${qdecproject_SRCS})

