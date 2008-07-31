project(mris_inflate)
include_directories(${mris_inflate_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_inflate_SRCS
mris_inflate.c
)


add_executable(mris_inflate ${mris_inflate_SRCS})
target_link_libraries(mris_inflate ${FS_LIBS})
install(TARGETS mris_inflate DESTINATION bin)	

