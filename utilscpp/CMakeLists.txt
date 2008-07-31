project(utilscpp)
include_directories(${utilscpp_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(utilscpp_SRCS
face.cpp
fastloop.cpp
globals.cpp
loop.cpp
mris_fastmarching.cpp
mris_topology.cpp
patchdisk.cpp
segment.cpp
surface.cpp
vertex.cpp
)


add_library(utilscpp ${utilscpp_SRCS})
add_subdirectory(test)
