project(repair_siemens_file)
include_directories(${repair_siemens_file_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(repair_siemens_file_SRCS
repair_siemens_file.c
)


add_executable(repair_siemens_file ${repair_siemens_file_SRCS})
target_link_libraries(repair_siemens_file ${FS_LIBS})
install(TARGETS repair_siemens_file DESTINATION bin)	

