project(dummy)
include_directories(${dummy_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(dummy_SRCS
dummy.c
myown.c
)


add_executable(dummy ${dummy_SRCS})
target_link_libraries(dummy ${FS_LIBS})
install(TARGETS dummy DESTINATION bin)	

