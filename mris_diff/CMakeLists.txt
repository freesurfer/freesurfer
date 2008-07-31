project(mris_diff)
include_directories(${mris_diff_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_diff_SRCS
mris_diff.c
)


add_executable(mris_diff ${mris_diff_SRCS})
target_link_libraries(mris_diff ${FS_LIBS})
install(TARGETS mris_diff DESTINATION bin)	

