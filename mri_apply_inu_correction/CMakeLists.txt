project(mri_apply_inu_correction)
include_directories(${mri_apply_inu_correction_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_apply_inu_correction_SRCS
mri_apply_inu_correction.c
)


add_executable(mri_apply_inu_correction ${mri_apply_inu_correction_SRCS})
target_link_libraries(mri_apply_inu_correction ${FS_LIBS})
install(TARGETS mri_apply_inu_correction DESTINATION bin)	

