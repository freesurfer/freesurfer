project(mris_niters2fwhm)
include_directories(${mris_niters2fwhm_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mris_niters2fwhm_SRCS
mris_niters2fwhm.c
)


add_executable(mris_niters2fwhm ${mris_niters2fwhm_SRCS})
target_link_libraries(mris_niters2fwhm ${FS_LIBS})
install(TARGETS mris_niters2fwhm DESTINATION bin)	

