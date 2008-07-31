project(mri_extract_ma_intensity)
include_directories(${mri_extract_ma_intensity_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_extract_ma_intensity_SRCS
mri_extract_ma_intensity.c
)


add_executable(mri_extract_ma_intensity ${mri_extract_ma_intensity_SRCS})
target_link_libraries(mri_extract_ma_intensity ${FS_LIBS})
install(TARGETS mri_extract_ma_intensity DESTINATION bin)	

