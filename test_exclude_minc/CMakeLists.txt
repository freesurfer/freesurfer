project(test_exclude_minc)
include_directories(${test_exclude_minc_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(test_exclude_minc_SRCS
test_exclude_minc.c
)


add_executable(test_exclude_minc ${test_exclude_minc_SRCS})
target_link_libraries(test_exclude_minc ${FS_LIBS})
install(TARGETS test_exclude_minc DESTINATION bin)	
