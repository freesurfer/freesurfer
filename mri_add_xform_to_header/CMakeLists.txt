project(mri_add_xform_to_header)
include_directories(${mri_add_xform_to_header_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_add_xform_to_header_SRCS
mri_add_xform_to_header.c
)


add_executable(mri_add_xform_to_header ${mri_add_xform_to_header_SRCS})
target_link_libraries(mri_add_xform_to_header ${FS_LIBS})
install(TARGETS mri_add_xform_to_header DESTINATION bin)	

