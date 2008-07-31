project(label2patch)
include_directories(${label2patch_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(label2patch_SRCS
label2patch.c
)


add_executable(label2patch ${label2patch_SRCS})
target_link_libraries(label2patch ${FS_LIBS})
install(TARGETS label2patch DESTINATION bin)	

