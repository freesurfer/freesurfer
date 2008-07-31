project(hiam_make_template)
include_directories(${hiam_make_template_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(hiam_make_template_SRCS
hiam_make_template.c
)


add_executable(hiam_make_template ${hiam_make_template_SRCS})
target_link_libraries(hiam_make_template ${FS_LIBS})
install(TARGETS hiam_make_template DESTINATION bin)	

