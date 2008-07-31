project(mri_make_template)
include_directories(${mri_make_template_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_make_template_SRCS
mri_make_template.c
)


add_executable(mri_make_template ${mri_make_template_SRCS})
target_link_libraries(mri_make_template ${FS_LIBS})
install(TARGETS mri_make_template DESTINATION bin)	

