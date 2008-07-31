project(mris_add_template)
include_directories(${mris_add_template_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_add_template_SRCS
mris_add_template.c
)


add_executable(mris_add_template ${mris_add_template_SRCS})
target_link_libraries(mris_add_template ${FS_LIBS})
install(TARGETS mris_add_template DESTINATION bin)	
