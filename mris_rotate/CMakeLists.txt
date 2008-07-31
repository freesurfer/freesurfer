project(mris_rotate)
include_directories(${mris_rotate_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_rotate_SRCS
mris_rotate.c
)


add_executable(mris_rotate ${mris_rotate_SRCS})
target_link_libraries(mris_rotate ${FS_LIBS})
install(TARGETS mris_rotate DESTINATION bin)	


