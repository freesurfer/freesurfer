project(mris_flatten)
include_directories(${mris_flatten_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_flatten_SRCS
mris_flatten.c
)


add_executable(mris_flatten ${mris_flatten_SRCS})
target_link_libraries(mris_flatten ${FS_LIBS})
install(TARGETS mris_flatten DESTINATION bin)	

