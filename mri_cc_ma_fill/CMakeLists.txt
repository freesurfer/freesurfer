project(mri_cc_ma_fill)
include_directories(${mri_cc_ma_fill_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_cc_ma_fill_SRCS
mri_cc_ma_fill.c
)


add_executable(mri_cc_ma_fill ${mri_cc_ma_fill_SRCS})
target_link_libraries(mri_cc_ma_fill ${FS_LIBS})
install(TARGETS mri_cc_ma_fill DESTINATION bin)	

