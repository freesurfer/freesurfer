project(tridec)
include_directories(${tridec_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(tridec_SRCS
tridec.c
)


add_executable(tridec ${tridec_SRCS})
target_link_libraries(tridec ${FS_LIBS})
install(TARGETS tridec DESTINATION bin)	

