project(mri_dct_align)
include_directories(${mri_dct_align_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(mri_dct_align_SRCS
mri_dct_align.c
mri_dct_align_binary.c
)


add_executable(mri_dct_align ${mri_dct_align_SRCS})
target_link_libraries(mri_dct_align ${FS_LIBS})
install(TARGETS mri_dct_align DESTINATION bin)	

