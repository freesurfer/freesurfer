project(mris_label2annot)
include_directories(${mris_label2annot_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_label2annot_SRCS
mris_label2annot.c
)


add_executable(mris_label2annot ${mris_label2annot_SRCS})
target_link_libraries(mris_label2annot ${FS_LIBS})
install(TARGETS mris_label2annot DESTINATION bin)	

