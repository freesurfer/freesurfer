project(mri_cc_medial_axis)
include_directories(${mri_cc_medial_axis_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_cc_medial_axis_SRCS
mri_cc_medial_axis.c
)


add_executable(mri_cc_medial_axis ${mri_cc_medial_axis_SRCS})
target_link_libraries(mri_cc_medial_axis ${FS_LIBS})
install(TARGETS mri_cc_medial_axis DESTINATION bin)	

