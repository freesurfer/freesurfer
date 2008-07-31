project(template)
include_directories(${template_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(template_SRCS
main_template.c
)


add_executable(template ${template_SRCS})
target_link_libraries(template ${FS_LIBS})
install(TARGETS template DESTINATION bin)	


