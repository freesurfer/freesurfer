project(mris_segment_vals)
include_directories(${mris_segment_vals_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_segment_vals_SRCS
mris_segment_vals.c
)


add_executable(mris_segment_vals ${mris_segment_vals_SRCS})
target_link_libraries(mris_segment_vals ${FS_LIBS})
install(TARGETS mris_segment_vals DESTINATION bin)	


