project(mris_translate_annotation)
include_directories(${mris_translate_annotation_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_translate_annotation_SRCS
mris_translate_annotation.c
)


add_executable(mris_translate_annotation ${mris_translate_annotation_SRCS})
target_link_libraries(mris_translate_annotation ${FS_LIBS})
install(TARGETS mris_translate_annotation DESTINATION bin)	


