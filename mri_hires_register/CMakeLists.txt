project(mri_hires_register)
include_directories(${mri_hires_register_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_hires_register_SRCS
mri_compute_distances.c
mri_hires_register.c
mri_linear_align.c
mri_linear_align_binary.c
mri_nl_align.c
mri_nl_align_binary.c
)


add_executable(mri_hires_register ${mri_hires_register_SRCS})
target_link_libraries(mri_hires_register ${FS_LIBS})
install(TARGETS mri_hires_register DESTINATION bin)	

