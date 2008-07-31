project(optseq2)
include_directories(${optseq2_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(optseq2_SRCS
optseq2.c
)


add_executable(optseq2 ${optseq2_SRCS})
target_link_libraries(optseq2 ${FS_LIBS})
install(TARGETS optseq2 DESTINATION bin)	


