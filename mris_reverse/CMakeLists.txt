project(mris_reverse)
include_directories(${mris_reverse_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_reverse_SRCS
mris_reverse.c
)


add_executable(mris_reverse ${mris_reverse_SRCS})
target_link_libraries(mris_reverse ${FS_LIBS})
install(TARGETS mris_reverse DESTINATION bin)	

