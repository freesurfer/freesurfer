project(x86cpucaps)
include_directories(${x86cpucaps_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(x86cpucaps_SRCS
check_sse.c
cpucaps.c
x86cpucaps.c
)


add_executable(x86cpucaps ${x86cpucaps_SRCS})
target_link_libraries(x86cpucaps ${FS_LIBS})
install(TARGETS x86cpucaps DESTINATION bin)	

