project(mris_fwhm)
include_directories(${mris_fwhm_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_fwhm_SRCS
mris_fwhm.c
)


add_executable(mris_fwhm ${mris_fwhm_SRCS})
target_link_libraries(mris_fwhm ${FS_LIBS})
install(TARGETS mris_fwhm DESTINATION bin)	

