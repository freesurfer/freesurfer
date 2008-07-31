project(rbftest)
include_directories(${rbftest_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(rbftest_SRCS
rbftest.c
)


add_executable(rbftest ${rbftest_SRCS})
target_link_libraries(rbftest ${FS_LIBS})
install(TARGETS rbftest DESTINATION bin)	


