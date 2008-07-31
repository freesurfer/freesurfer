project(mri_fwhm)
include_directories(${mri_fwhm_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_fwhm_SRCS
mri_fwhm.c
)


add_executable(mri_fwhm ${mri_fwhm_SRCS})
target_link_libraries(mri_fwhm ${FS_LIBS})
install(TARGETS mri_fwhm DESTINATION bin)	
