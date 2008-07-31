project(check_siemens_dir)
include_directories(${check_siemens_dir_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(check_siemens_dir_SRCS
check_siemens_dir.c
)


add_executable(check_siemens_dir ${check_siemens_dir_SRCS})
target_link_libraries(check_siemens_dir ${FS_LIBS})
install(TARGETS check_siemens_dir DESTINATION bin)	

