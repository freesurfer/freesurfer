project(histo_register_block)
include_directories(${histo_register_block_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(histo_register_block_SRCS
histo_register_block.c
)


add_executable(histo_register_block ${histo_register_block_SRCS})
target_link_libraries(histo_register_block ${FS_LIBS})
install(TARGETS histo_register_block DESTINATION bin)	

