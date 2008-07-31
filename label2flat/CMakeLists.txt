project(label2flat)
include_directories(${label2flat_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(label2flat_SRCS
label2flat.c
)


add_executable(label2flat ${label2flat_SRCS})
target_link_libraries(label2flat ${FS_LIBS})
install(TARGETS label2flat DESTINATION bin)	
