project(mris_make_template)
include_directories(${mris_make_template_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_make_template_SRCS
mris_make_template.c
)


add_executable(mris_make_template ${mris_make_template_SRCS})
target_link_libraries(mris_make_template ${FS_LIBS})
install(TARGETS mris_make_template DESTINATION bin)	

