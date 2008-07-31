project(mris_expand)
include_directories(${mris_expand_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_expand_SRCS
mris_expand.c
)


add_executable(mris_expand ${mris_expand_SRCS})
target_link_libraries(mris_expand ${FS_LIBS})
install(TARGETS mris_expand DESTINATION bin)	
