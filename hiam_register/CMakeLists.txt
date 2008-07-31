project(hiam_register)
include_directories(${hiam_register_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(hiam_register_SRCS
hiam_register.c
)


add_executable(hiam_register ${hiam_register_SRCS})
target_link_libraries(hiam_register ${FS_LIBS})
install(TARGETS hiam_register DESTINATION bin)	
