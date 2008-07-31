project(tkregister2)
include_directories(${tkregister2_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(tkregister2_SRCS
tkregister2.c
)


add_executable(tkregister2 ${tkregister2_SRCS})
target_link_libraries(tkregister2 ${FS_LIBS})
install(TARGETS tkregister2 DESTINATION bin)	

add_subdirectory(cmdl)
