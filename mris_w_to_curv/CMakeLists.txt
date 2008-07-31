project(mris_w_to_curv)
include_directories(${mris_w_to_curv_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_w_to_curv_SRCS
mris_w_to_curv.c
)


add_executable(mris_w_to_curv ${mris_w_to_curv_SRCS})
target_link_libraries(mris_w_to_curv ${FS_LIBS})
install(TARGETS mris_w_to_curv DESTINATION bin)	


