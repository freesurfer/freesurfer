project(mris_seg2annot)
include_directories(${mris_seg2annot_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_seg2annot_SRCS
mris_seg2annot.c
)


add_executable(mris_seg2annot ${mris_seg2annot_SRCS})
target_link_libraries(mris_seg2annot ${FS_LIBS})
install(TARGETS mris_seg2annot DESTINATION bin)	

